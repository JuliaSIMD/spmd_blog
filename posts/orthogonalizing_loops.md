+++
title = "Orthogonalize Indexes"
hascode = true
date = Date(2022, 3, 14)
rss = "This post goes over the significance of unimodular matrices in the
context of loop optimizations"
+++

@def tags = ["unimodular", "loop-transformations"]

In a convolutional neural network, we are likely to have the following two
pieces of code, calculating the forward pass, and then the revese pass for
calculating the gradient with respect to `A`.
```julia-repl
julia> using OffsetArrays

julia> function conv!(_C, _A, _B)
    A = OffsetArray(_A, OffsetArrays.Origin(0,0))
    B = OffsetArray(_B, OffsetArrays.Origin(0,0))
    C = OffsetArray(_C, OffsetArrays.Origin(0,0))
    I, J = size(B)
    M, N = size(C)
    for n = 0:N-1, m = 0:M-1, i = 0:I-1, j = 0:J-1
        C[m,n] += A[m + i, n + j] * B[i, j]
    end
    return C
end

julia> function convpullback!(_Ā, _B, _C̄)
    Ā = OffsetArray(_Ā, OffsetArrays.Origin(0,0))
    B = OffsetArray(_B, OffsetArrays.Origin(0,0))
    C̄ = OffsetArray(_C̄, OffsetArrays.Origin(0,0))
    I, J = size(B)
    M, N = size(C̄)
    for n = 0:N-1, m = 0:M-1, i = 0:I-1, j = 0:J-1
        Ā[m + i, n + j] += B[i, j] * C̄[m, n]
    end
    return Ā
end
```
We left off the code for the reverse pass with respect to `B`, as it doesn't
feature the problem this blog post focuses on solving: 
```julia
Ā[m + i, n + j] += 
```

While the forward pass is easy to optimize via register tiling, this is not
the case for the pullback: the index into $\bar{A}$ is dependent on all four loop
induction variables, making it impossible to hoist these loads and stores
out of any loops, forcing us to reload and restore memory on every iteration,
meaning it takes many times more CPU instructions to evaluate.

I'll add a post detailing register tiling, but for now just note that asside from
reducing the loads and stores to $\bar{A}$ by factors equal to the loops they're 
hoisted out of, it would also results in several times fewer loads from `B` and
from $\bar{C}$, enabling code to start attaining significant fractions of peak flops.
An order of magnitude difference in performance is often a fair ballpark for
the benefit of register tiling.

So, now the question that's the focus of this post: can we re-index the memory
accesses so that $\bar{A}$ is dependent on only two loops?
That is, we want to produce a new set of loops that looks more like
```julia
for w in W, x in X
    Āwx = Ā[w, x]
    for y in Y, z in Z
        Āwx += B[???] * C̄[???]
    end
    Ā[w, x] = Āwx
end
```
which would allow us to register tile. If we can, what are the new ranges
`W`, `X`, `Y`, and `Z`, and what are the new indices into `B` and $\bar{C}$?
Can we develop a general algorithm?

Often, a human can look at a problem and reason their way to a solution without
too much difficulty. While our goal is to create a general algorithm to remove
the need for expert human intervention (especially to enable codegen tools like
reverse diff on loops to be produce optimal code; it is an express goal that
naive loops + an AD tool like [Enzyme](https://github.com/EnzymeAD/Enzyme) + the new LoopVectorization
will achieve or best state of the art performance across a wide range of loops),
I find this is a useful starting point for building an intuition of the problem,
which in turn can help point to general algorithms.

We want to set `w = m + i`, and `x = n + j`, thus `W` must iterate over the full
range of values attained by `m + i`, i.e. `w = 0:M+N-2`, and `x = 0:N+J-2`.
We may now naively try setting the indices of `B[i,j]` to `B[y,z]` and working
through the implications; what would this imply about the indicesof $\bar{C}$, and
about the ranges `Y` and `Z`?

First, it's straightforward to find that we must have `C̄[w-y, x-z]`, as
`m = w - i = w - y`, and we can prove `n` similarly.

For the bounds on `y`, note that `y = i`, and `0 <= i <= I-1`, thus
`0 <= y <= I-1`. However, we also have that `0 <= m <= M-1`, therefore
`0 <= w-y <= M-1`, which means `y <= w` and `y >= w - (M-1)`. Taking
the intersection yields `max(0, w - (M-1)) <= y <= min(I-1, w)`.
We can apply the same argument for `z` to produce the bounds
`max(0, x - (N-1)) <= z <= min(J-1, x)`.

To make this algorithmic, we represent the loop bounds via a system of
inequalities:
```julia
[ 1  0  0  0    [ m         [ M-1
 -1  0  0  0      n    .<=     0
  0  1  0  0      i           N-1
  0 -1  0  0      j ]          0
  0  0  1  0                  I-1
  0  0 -1  0                   0
  0  0  0  1                  J-1
  0  0  0 -1 ]                 0 ]
# A * [m;n;i;j] <= b
```
and the indices with a matrix
```julia
[ 1 0 1 0   * [ m         [ m + i
  0 1 0 1       n    .==    n + j
  1 0 0 0       i           m
  0 1 0 0       j ]         n
  0 0 1 0                   i
  0 0 0 1 ]                 j     ]

# X * [m;n;i;j] == [index expressions...]
```
Letting the index matrix be `X`, we can then frame our objective as
finding some invertible matrix `K` such that the first two rows of `X*K`
are linearly independent single-element vectors, with `1` as the single element.
We can permute the columns of `K` arbitrarily to simplify this to say `X*K`'s first
two rows should be `hcat(I(2), 0)`.
With this, we can insert $\textbf{I} = \textbf{KK}^{-1}$ to apply the transform.
This also defines our new loop induction variables
$
\textbf{K}^{-1}\begin{bmatrix}m\\n\\i\\j\end{bmatrix}=
\begin{bmatrix}w\\x\\y\\z\end{bmatrix}
$.
Define $\textbf{w}$ as our new vector of induction variables, we can also express
our new loop inequalities simply as $\textbf{AK}^{-1}\textbf{v} <= \textbf{b}$.

(Note that $\textbf{K}^{-1}$ is a simple representation of a linear loop schedule
-- it is a linear function that gives the order in which iterations of the loops
are evaluated. We will discuss affine schedules in much greater detail in a future
post.)


This means that all we have left to do is actually find the matrix $K$.
As $K$ must be both an integer matrix and invertible, it is unimodular.
As the first two rows of $\textbf{XK}$ correspond to the identity matrix,
we know the first two rows of $X$ equal the first two rows of $\textbf{K}^{-1}$.

From here, we realize that we can find such a matrix by modifying the algorithms
typically used to bring matrices into reduced forms, such as reduced echelon form
or the Hermite normal form. We can use column pivots and row additions (using 
integer multipliers), as these preserve the determinant and maintain the status
as an integer matrix).
We diagonalize one row at a time. When we encounter a row we cannot diagonalize, that 
means this index cannot produce a unimodular matrix in combination with the preceding
indices, thus we reject it and move on to the next index. In this way our matrix $\textbf{K}^{-1}$
can be a subset of up to `size(X,2)` rows of $\textbf{X}$, not just the first two.
The important characteristics here are that it allows us to choose which indices to
prioritize -- e.g. those that would enable register tiling -- while still simplifying/
maintaining simplicity in some other indices.
If we run out of rows of $\textbf{X}$, then the algorithm will still succeed, but then
$\textbf{K}^{-1}$ will include rows not present in $\textbf{X}$.

Inputting the loop bounds (as `0 <= i_0 <= M-1`, `0 <= i_1 <= N-1`, `0 <= i_2 <= O-1`,
and `0 <= i_3 <= P-1`) and the indices, we get the following output post-transformation:
```
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <=  ( M + O - 2 )
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <=  ( N + P - 2 )
Loop 2 lower bounds:
i_2 >=  ( - M + 1 )  + i_0
i_2 >= 0
Loop 2 upper bounds:
i_2 <= i_0
i_2 <=  ( O - 1 )
Loop 3 lower bounds:
i_3 >=  ( - N + 1 )  + i_1
i_3 >= 0
Loop 3 upper bounds:
i_3 <= i_1
i_3 <=  ( P - 1 )

New ArrayReferences:
ArrayReference 0 (dim = 2):
{ Induction Variable: 0 }
 ( M + O - 1 )  * ({ Induction Variable: 1 })



ArrayReference 1 (dim = 2):
{ Induction Variable: 2 }
 ( O )  * ({ Induction Variable: 3 })



ArrayReference 2 (dim = 2):
{ Induction Variable: 0 } - { Induction Variable: 2 }
 ( M )  * ({ Induction Variable: 1 } - { Induction Variable: 3 })
```
Our desired/expected outcome. We can now register tile the loop.



Now let us move to a more artificial example, a deliberately badly written matrix-multiply:
```julia-repl
julia> using OffsetArrays

julia> function badmul!(C,A,B)
           M,N = size(C)
           K = size(B,1); fill!(C,0)
           for i in 0:M+N+K-3, l in max(0,i+1-N):min(M+K-2,i), j in max(0,l+1-K):min(M-1,l)
               C[j,i-l] += A[j,l-j]*B[l-j,i-l]
           end
           C
       end
badmul! (generic function with 1 method)

julia> M,K,N = 5,6,7
(5, 6, 7)

julia> A = OffsetArray(rand(M,K),-1,-1);

julia> B = OffsetArray(rand(K,N),-1,-1);

julia> C = OffsetArray(rand(M,N),-1,-1);

julia> badmul!(C,A,B)
5×7 OffsetArray(::Matrix{Float64}, 0:4, 0:6) with eltype Float64 with indices 0:4×0:6:
 1.9536   1.20026   2.20549  1.11528   1.77055  2.14641  1.75909
 1.54733  1.10697   1.6769   1.13867   1.06312  1.71708  1.65325
 1.76249  0.908232  1.63373  0.995246  1.19762  1.69865  1.71706
 1.46832  0.939424  1.78488  1.36437   1.09265  1.69105  1.35382
 1.27257  0.770909  1.14775  0.596495  0.92855  1.21593  1.26251

julia> parent(A)*parent(B)
5×7 Matrix{Float64}:
 1.9536   1.20026   2.20549  1.11528   1.77055  2.14641  1.75909
 1.54733  1.10697   1.6769   1.13867   1.06312  1.71708  1.65325
 1.76249  0.908232  1.63373  0.995246  1.19762  1.69865  1.71706
 1.46832  0.939424  1.78488  1.36437   1.09265  1.69105  1.35382
 1.27257  0.770909  1.14775  0.596495  0.92855  1.21593  1.26251
```
Entering the `badmul!` loopnest into the orthoganlizer produces:
```
Skewed loop nest:
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <=  ( M - 1 )
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <=  ( N - 1 )
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <=  ( O - 1 )

New ArrayReferences:
ArrayReference 0 (dim = 2):
{ Induction Variable: 0 }
 ( M )  * ({ Induction Variable: 1 })



ArrayReference 1 (dim = 2):
{ Induction Variable: 0 }
 ( O )  * ({ Induction Variable: 2 })



ArrayReference 2 (dim = 2):
{ Induction Variable: 2 }
 ( M )  * ({ Induction Variable: 1 })
```
Recovering a perfectly-normal tripple-loop matrix multiply, that is straightforward to
optimize via register and cache tiling.

<!---
Let `L` be the vector of loop ind vars, so `L = [i, l, j]`, and `I_{X}` be the
set of indices for array `{X}`. Then let `I_{X} = A_{X} * L`
```julia-repl
julia> A_A = [ 0  0  1;
               0  1 -1];

julia> A_B = [ 0  1 -1;
               1 -1  0];

julia> A_C = [ 0  0  1;
               1 -1  0];
```
What we want is that each row has only a single non-zero element.
Something something basis.

```julia-repl
julia> U = reduce(hcat, unique(Iterators.flatten(Iterators.map(eachrow, (A_A, A_B, A_C)))))'
3×3 adjoint(::Matrix{Int64}) with eltype Int64:
 0   0   1
 0   1  -1
 1  -1   0
 
julia> lu(Rational.(U))
LU{Rational{Int64}, Matrix{Rational{Int64}}, Vector{Int64}}
L factor:
3×3 Matrix{Rational{Int64}}:
 1//1  0//1  0//1
 0//1  1//1  0//1
 0//1  0//1  1//1
U factor:
3×3 Matrix{Rational{Int64}}:
 1//1  -1//1   0//1
 0//1   1//1  -1//1
 0//1   0//1   1//1

julia> inv(ans)
3×3 Matrix{Rational{Int64}}:
 1//1  1//1  1//1
 1//1  1//1  0//1
 1//1  0//1  0//1
```
We have `L = A_{X} \ I_{X}`, but `A_{X}` is of course not necessarily square.
Let `U*L = O` so that `L = U^{-1}*O`
What we really want to do is transform the loops. Loop bounds are

\begin{align}
\textbf{W}*\textbf{L} &\le \textbf{B}\\
(\textbf{W}*\textbf{U}^{-1})*\textbf{O} &\le \textbf{B}
\end{align}

In our example, we have

\begin{align}
0 &\le i \le M+N+K-3\\
0 &\le l \le M+K-2\\
i-(N-1) &\le l \le i\\
0 &\le j \le M-1\\
l-(K-1) &\le j \le l
\end{align}

```julia-repl
julia> W = [ 1  0  0    #  [ i       [ M+N+K-3
            -1  0  0    #    l    <=   0
             0  1  0    #    j ]       M+K-2
             0 -1  0    #              0
             1 -1  0    #              N-1
            -1  1  0    #              0
             0  0  1    #              M-1
             0  0 -1    #              0
             0  1 -1    #              K-1
             0 -1  1 ]  #              0      ]

julia> WAi = Int.(W / lu(Rational.(U)))
10×3 Matrix{Int64}: # let O = [x, y, z]
  1   1   1  #  [ x       [ M+N+K-3
 -1  -1  -1  #    y    <=   0
  1   1   0  #    z ]       M+K-2
 -1  -1   0  #              0
  0   0   1  #              N-1
  0   0  -1  #              0
  1   0   0  #              M-1
 -1   0   0  #              0
  0   1   0  #              K-1
  0  -1   0  #              0      ]
```
So this tells us:

\begin{align}
0 &\le z \le N-1\\
0 &\le x \le M-1\\
0 &\le y \le K-1
\end{align}

Now, what happens if we have loops more like
```julia
for i in 0:I-1, j in 0:J-1, k in 0:K-1, l in 0:L-1
    C[i,j] += A[i+k,j+l]*B[k,l]
end
```
more interesting, then...
Lets say someone runs reverse diff on this
With respect to the kernel:
```julia
for i in 0:I-1, j in 0:J-1, k in 0:K-1, l in 0:L-1
  ̄B[k,l] += A[i+k,j+l]*̄C[i,j]
end
```
This is nice and easy to optimize. Let's move on.
With respect to whatever was being convolved with the kernel:
```julia
for i in 0:I-1, j in 0:J-1, k in 0:K-1, l in 0:L-1
    Ā[i+k,j+l] += C̄[i,j]*B[k,l]
end
```
(You may need this if you have a convolution like this that isn't the first
layer of a conv net, as it must backprop $\bar{A}$ to the preceding layer)
Let's try and get it so we can actually hoist $\bar{A}$ out of a few loops and
tile.

```julia
A_̄A = [1 0 1 0;
       0 1 0 1];

A_̄C = [1 0 0 0;
       0 1 0 0];

A_B    = [0 0 1 0;
          0 0 0 1];

```

```julia
for ik in 0:I+K-2, jl in 0:J+L-2, k in max(0,ik-(I+K-2)):min(K-1,ik), l in max(0,jl-(J+L-2)):min(L-1,jl)
  ̄A[ik, jl] += ̄C[ik-k,jl-l]*B[k,l]
end
U = [ 1 0 1 0    # [ i   = [ ik
      0 1 0 1    #   j       jl
      0 0 1 0    #   k       k
      0 0 0 1 ]  #   l ]     l ]
# L = [i, j, k, l]
# N = [ik, jl, k, l]
#
# W * L <= B
# U * L =  N
# L = U^{-1}*N
# (W*U^{-1}) * N <= B
julia> Uinv = Int.(inv(lu(Rational.(U))))
4×4 Matrix{Int64}:
 1  0  -1   0
 0  1   0  -1
 0  0   1   0
 0  0   0   1

# 0 <= i <= I-1
# 0 <= j <= J-1
# 0 <= k <= K-1
# 0 <= l <= L-1
W = [ 1  0  0  0    #  [ i    <= [ I-1
     -1  0  0  0    #    j         0
      0  1  0  0    #    k         J-1
      0 -1  0  0    #    l ]       0
      0  0  1  0    #              K-1
      0  0 -1  0    #              0
      0  0  0  1    #              L-1
      0  0  0 -1 ]  #              0   ]

julia> W*Uinv
8×4 Matrix{Int64}:
  1   0  -1   0    #  [ ik   <= [ I-1
 -1   0   1   0    #    jl        0
  0   1   0  -1    #    k         J-1
  0  -1   0   1    #    l ]       0
  0   0   1   0    #              K-1
  0   0  -1   0    #              0
  0   0   0   1    #              L-1
  0   0   0  -1 ]  #              0   ]

# 0 <= k <= K-1
# 0 <= l <= L-1
# k <= ik <= I-1 + k
# l <= jl <= J-1 + l
# now, we want to permute the `ik` loop to the outside
# 0 <= ik <= I-1 + K-1 = I+K-2
# 0 <= k <= K-1
# ik - I-1 <= k <= ik
# now, we want to permute the `jl` loop to the outside
# 0 <= jl <= J-1 + L-1 = J+L-2
# 0 <= l <= L-1
# jl-(J-1) <= l <= jl
```

We can take the approach of starting with
\begin{align}
\begin{bmatrix}
1 & 0 & 1 & 0\\
0 & 1 & 0 & 1\\
\end{bmatrix}
\begin{bmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
\end{bmatrix}
\end{align}
try to reduce, and then keep adding rows from other indices one at a time until full rank is achieved.
If an added row is a linear combination of existing rows, drop it and revert changes.



We should solve for possible solutions.
Reducing `U'` to reduced row echelon form, we get

\begin{align}
\begin{bmatrix}
1 & 0 & 1 & 0 & 0 & 0\\
0 & 1 & 0 & 1 & 0 & 0\\
1 & 0 & 0 & 0 & 1 & 0\\
0 & 1 & 0 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
\rightarrow\\
\begin{bmatrix}
1 & 0 & 0 & 0 & 1 & 0\\
0 & 1 & 0 & 0 & 0 & 1\\
0 & 0 & 1 & 0 & -1 & 0\\
0 & 0 & 0 & 1 & 0 & -1
\end{bmatrix}
\begin{bmatrix}
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1\\
1 & 0 & -1 & 0\\
0 & 1 & 0 & -1
\end{bmatrix}
\end{align}
Letting `O = [w, x, y, z]` be our new indices, we can interpret this as meaning we must have
\begin{align}


\end{align}



Given `a`, `b`, `c`, and solve for `x`, `y`, `z`
[a b c]
[x y z]

=>

[x y z]            = [1 y z]
[kx+a, ky+b, kz+c] = [0 i j] such that `gcd(i, j) == 1`
x == 1 => k = -a
gcd(-ay+b, -az+c) = 1
gcd(b, -az+c) = 1

-->


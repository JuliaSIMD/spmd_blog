+++
title = "Unimodular Matrices and Loops"
hascode = true
date = Date(2022, 3, 2)
rss = "This post goes over the significance of unimodular matrices in the
context of loop optimizations"
+++

@def tags = ["unimodular", "loop-transformation"]

A unimodular matrix is a square integer matrix with determinate $\pm 1$. With
some basic linear algebra facts, we can show that a square integer matrix is
unimodular if and only if its inverse exists and the inverse matrix is also an
integer matrix.

Let $A$ be an unimodular matrix. As the determinate of an integer matrix must
also be an integer, $\adj(A)$ is also an integer matrix. It follows that $A^{-1}
= \frac{1}{\det{A}} \adj(A)$ is an integer matrix.

To show the converse, let $A$ be an integer square matrix such that $A^{-1}$ is
an integer matrix as well. Since $\det(A) \det(A^{-1}) = \det(A A^{-1}) = \det(I)
= 1$ and $\det(A), \det(A^{-1}) \in \ZZ$, we can conclude that $A$ and $A^{-1}$
are unimodular.

Say we have code like the following:
```julia-repl
julia> using OffsetArrays

julia> function badmul2!(C,A,B)
           M,N = size(C)
           K = size(B,1); fill!(C,0)
           for i in 0:M+N+K-3, l in max(0,i+1-N):min(M+K-2,i), j in max(0,l+1-K):min(M-1,l)
               C[j,i-l] += A[j,l-j]*B[l-j,i-l]
           end
           C
       end
badmul2! (generic function with 1 method)

julia> M,K,N = 5,6,7
(5, 6, 7)

julia> A = OffsetArray(rand(M,K),-1,-1);

julia> B = OffsetArray(rand(K,N),-1,-1);

julia> C = OffsetArray(rand(M,N),-1,-1);

julia> badmul2!(C,A,B)
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
    Bbar[k,l] += A[i+k,j+l]*Cbar[i,j]
end
```
This is nice and easy to optimize. Let's move on.
With respect to whatever was being convolved with the kernel:
```julia
for i in 0:I-1, j in 0:J-1, k in 0:K-1, l in 0:L-1
    Abar[i+k,j+l] += Cbar[i,j]*B[k,l]
end
```
(You may need this if you have a convolution like this that isn't the first
layer of a conv net, as it must backprop Abar to the preceding layer)
Let's try and get it so we can actually hoist `Abar` out of a few loops and
tile.

```julia
A_Abar = [1 0 1 0;
          0 1 0 1];

A_Cbar = [1 0 0 0;
          0 1 0 0];

A_B    = [0 0 1 0;
          0 0 0 1];

```

```julia
for ik in 0:I+K-2, jl in 0:J+L-2, k in max(0,ik-(I+K-2)):min(K-1,ik), l in max(0,jl-(J+L-2)):min(L-1,jl)
    Abar[ik, jl] += Cbar[ik-k,jl-l]*B[k,l]
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









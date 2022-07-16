+++
title = "Dependence Polyhedra"
hascode = true
date = Date(2022, 3, 28)
rss = "polyhedra dependency"
+++

@def tags = ["dependence polyhedra", "polyhedral optimization", "loop-transformations"]

Dependency Polyhedra

Ready for some math dump without explanation or context?
All right, let's go!!

Any errors can be assumed to be due to a pedagogical focus rather than lack of
rigor or a mistake on our part. But please file an issue or PR to let us know
anyway.

We have this loop:
```julia
for i = 1:N, j = 1:i
    A[i,j] = foo(A[i,i])
end
```

$
1 \le i_0 \le N\\
1 \le j_0 \le i_0\\
1 \le i_1 \le N\\
1 \le j_1 \le i_1\\
i_0 = i_1\\
j_0 = i_1
$

which simplifies to

$
i_0 \le N\\
1 \le j_1\\
j_0 \le i_0\\
j_1 \le i_1\\
i_0 = i_1\\
j_1 = i_1
$

Given two memory references, you get the polyhedra for the iteration spaces of
both, which gives you a higher dimensional space, and then you have equalities
to match the indices that lead to identical memory addresses in the accessed
arrays. These can restrict the space.

It's straightforward/algorithmic to derive this sort of thing from a program
(de-linearization of indices aside). However, when we have two accesses that are
not both reads to the same memory address, we must also know which of the two
came first. If we have two reads, this may still be of interest, as minimizing
the distance between both reads can increase cache reuse, but switching their
order does not change correctness/the program's results, so we don't add any
constraints based on pairs of reads.

How do we figure out which one came first?

Well, sometimes neither is first, because their isn't actually a dependence at
all. E.g.
```julia
for i = 1:N, j = 1:i-1
    A[i,j] = foo(A[j,i])
end
```
the associated dependence polyhedra is empty. Proof is left as an exercise for
the reader. Note this is one of LoopModel's test cases, so if you get stuck,
[check out our C++ code](https://github.com/JuliaSIMD/LoopModels/blob/e779fc966b37fb71c1880366310e491734f75118/test/dependence_test.cpp#L195)
as it's pristine readability will elucidate all.

But in non-empty cases, we use Farkas' lemma. There might be more
straightforward approaches to getting the direction, but we'll need to apply it
anyway, so making use of it for this proof adds little marginal computational
cost.

So what is Farkas' lemma? Read the other blog posts.

So, to find out whether the load or the store came first, we use Farkas lemma to
produce contradictory constraints, and then find out which one of these
constraints is satisfied by the original schedule of the code as written.

By telling Farkas that $0$ (i.e., the store) happens before $1$ (i.e., the load),
then Farkas gives us these two inequalities:

$
-\alpha_{i_0} - \alpha_{j_0} + \alpha_{i_1} + \alpha_{j_1} \le 0\\
-\alpha_{i_1} - \alpha_{j_0} + \alpha_{i_1} \le 0\\
$

Then, telling Farkas the opposite, it instead returns these two inequalities:

$
\alpha_{i_0} + \alpha_{j_0} - \alpha_{i_1} - \alpha_{j_1} \le 0\\
\alpha_{i_0} + \alpha_{j_0} - \alpha_{i_1} \le 0
$

So, how do we check which is violated?

These are constraints on valid schedules. The schedule of the code as written is
the identity plus some offsets that aren't important right now.
Read the blog post about schedules that will be written later as my schedule allows
for more information.

The important bit here is that we iterate over rows of this identity matrix
(intermixed with those constants indicating loop orders of imperfect nests).
For the first row, we have $\alpha_{i_0} = \alpha_{i_1} = 1$ and $\alpha_{j_0} =
\alpha_{j_1} = 0$.

Plugging that in, we get for our first equation:

$
-1 - 0 + 1 + 0 = 0 \le 0\\
-1 - 0 + 1 = 0 \le 0
$

and for our second

$
1 + 0 - 1 - 0 = 0 \le 0\\
1 + 0 - 1 = 0 \le 0
$

Now, for the second row of the identity matrix:

$
-0 -1 + 0 + 1 = 0 \le 0\\
-0 - 1 + 0 = -1 \le 0
$

and

$
0 + 1 - 0 - 1 = 0 \le 0\\
0 + 1 - 0 = 1 \le 0
$

Meaning this last constraint is invalid, and we have thus proved that the load
precedes the store.


Now, let's consider repeated memory accesses across time.
```julia
for i = 0:I-1, j = 0:J-1, k = 0:K-1
    A[i+j, j+k, i-k] = foo(A[i+j, j+k, i-k])
end
```
We can represent the mapping of loop induction variables to array memory
accesses using linear algebra by specifying a matrix $\textbf{R}$ where

$
\begin{bmatrix}1&1&0\\0&1&1\end{bmatrix} \begin{bmatrix}i\\j\\k\end{bmatrix}\\
= \begin{bmatrix} i+j\\j+k\end{bmatrix}
$

ModelingToolkit is build to help with exactly these sorts of problems (little
known fact: it also happens to have a bit of differential equation handling code)
```julia
julia> using ModelingToolkit

julia> A = [1  1  0
            0  1  1
			1  0 -1];

julia> ModelingToolkit.nullspace(A)
3×1 Matrix{Int64}:
 -1
  1
 -1
```
so, `[-i, j, -k]` is a "time direction" representing repeated accesses to the
same memory.
Now, we wonder, when holding the memory address constant, i.e. holding both $i+j$ and $j+k$
constant, in what direction does $-i + j -k$ change?
We need to know that to place constraints indicating whether one time is before
or after the other.

We can solve this via applying the same approach as above. We first get the direction of the dependencies when time is equal across both. Then, we pick one of the two directions, and see if the constraints assuming the reverse direction are satisfied. If they are, the direction was correct, else the time vector goes in the opposite direction.

Equality constraint simplification can remove any clear correspondence between constraints and time dimensions, thus results are not necessarilly easily interpretable. But the following list of citations suggests readers love lots of mathy output regardless of interpretability or context, so I'm presenting results anyway:

For the forward direction, the dependecy polyhedra:
$
-i_l - 2k_l + i_s + 2k_s = 0
-j_l - k_l + j_s + k_s = 0
-3k_l + 3k_s = 0
$

Schedule Constraints:
-v_1 + v_4 <= 0
-v_2 + v_5 <= 0
-v_0 + v_3 <= 0
-v_6 == 0

Bounding Constraints:
-v_7 <= 0
-v_8 <= 0
-v_9 <= 0
-v_10 <= 0
v_2 - v_5 - v_10 <= 0
v_0 - v_3 - v_8 <= 0
v_1 - v_4 - v_9 <= 0
-v_6 == 0


Reverse:
Dependence Poly x -> y:
v_0 <=  ( M - 1 )
-v_0 <= 0
v_1 <=  ( N - 1 )
-v_1 <= 0
v_2 <=  ( O - 1 )
-v_2 <= 0
v_3 <=  ( M - 1 )
-v_3 <= 0
v_4 <=  ( N - 1 )
-v_4 <= 0
v_5 <=  ( O - 1 )
-v_5 <= 0
-v_0 - 2v_2 + v_3 + 2v_5 == -1
-v_1 - v_2 + v_4 + v_5 == 0
-3v_2 + 3v_5 == -1

Schedule Constraints:
v_2 - v_5 <= 0
v_0 - v_3 <= 0
v_1 - v_4 <= 0
-v_6 == 0

Bounding Constraints:
-v_7 <= 0
-v_8 <= 0
-v_9 <= 0
-v_10 <= 0
-v_0 + v_3 - v_8 <= 0
-v_1 + v_4 - v_9 <= 0
-v_2 + v_5 - v_10 <= 0
-v_6 == 0
$





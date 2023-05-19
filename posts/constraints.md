+++
title = "Analyzing constraints"
hascode = true
date = Date(2023, 5, 17)
rss = "Analyzing loop constraints"
+++

@def tags = ["dependency analysis", "constraints", "loops"]

Let's say we're given a loop nest like:
```julia
using Test
function trisolve!(A, B, U)
  M, N = size(A)
  @assert size(B) == (M, N)
  @assert size(U) == (N, N)
  for m = 1:M
    for n = 1:N
      A[m, n] = B[m, n]
    end
    for n = 1:N
      A[m, n] /= U[n, n]
      for k = n+1:N
        A[m, k] -= A[m, n] * U[n, k]
      end
    end
  end
end
M, N = 4, 7
B = rand(M, N);
U = UpperTriangular(rand(N, N));
A = similar(B);
trisolve!(A, B, U)
@test A ≈ B / U # pass
```
When we analyze it in LLVM, scalar evolution presents us address evolution as a function of loop trip counts, that will essentially transform the loop into the one below:
```julia
using OffsetArrays, Test
function trisolve_preprocess!(_A, _B, _U)
  M, N = size(_A)
  @assert size(_B) == (M, N)
  @assert size(_U) == (N, N)
  A = OffsetArray(_A, -1, -1) # arrays all become 0-indexed
  B = OffsetArray(_B, -1, -1)
  U = OffsetArray(_U, -1, -1)
  for m = 0:M-1
    for n = 0:N-1
      A[m, n] = B[m, n]
    end
    for n = 0:N-1
      A[m, n] /= U[n, n]
      for k = 0:N-n-2 # loops all start at 0
        A[m, k+n+1] -= A[m, n] * U[n, k+n+1]
      end
    end
  end
end
M, N = 4, 7
B = rand(M, N);
U = UpperTriangular(rand(N, N));
A = similar(B);
trisolve_preprocess!(A, B, U)
@test A ≈ B / U # pass
```
That is, all arrays become 0-indexed (which is most convenient), and all loops start at 0 -- which is perhaps a reasonable canonicalization, but can cause us some difficulty when building up constraint systems.

Intuitively, both versions of the code represent the exact same program with the same traversal order of memory addresses, and thus should be interpreted the same way.

To understand an issue with the second representation let us first consider the dependence between the store `A[m,n] = A[m,n] / U[n,n]`, and the loads used for updating in the inner loop, that is `A[m, k]` for the first version:
```julia
A[m, k] = A[m, k] - A[m, n] * U[n, k]
```
or `A[m, k+n+1]` for the second version:
```julia
A[m, k+n+1] = A[m, k+n+1] - A[m, n] * U[n, k+n+1]
```
Let's take for granted that the loads at a particular address (`A[m, k]` and `A[m, k+n+1]`, respectively) must happen before the store `A[m,n]`. To see that this is the case, consider that once we compute `A[m, n]` for a particular value of `n`, say `n = x`, we update all values of `A[m, y]`, where `x < y` and `y` is still within bounds. Only on a later iteration after `n` has been incremented `y-x` times do we finally perform the `A[m,n]` store to that address. That is, the store happens after the load, because for any particular address, we stored after the load.

For the first example, we have the iteration spaces for `A[m,n]`, hereafter called the `target`, or `t`, and the load `A[m,k]`, henceforth referred to as the source, or `s`:
\begin{align}
\begin{bmatrix}
-1 & 0 & 0 & 1 & 0\\
-1 & 0 & 0 & 0 & 1\\
0 & 1 & 0 & -1 & 0\\
0 & 0 & 1 & 0 & -1\\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_t\\n_t
\end{bmatrix}\ge\textbf{0}\\
\begin{bmatrix}
-1 & 0 & 0 & 1 & 0 & 0\\
-1 & 0 & 0 & 0 & 1 & 0\\
-1 & 0 & 0 & 0 & -1 & 1\\
0 & 1 & 0 & -1 & 0 & 0\\
0 & 0 & 1 & 0 & -1 & 0\\
0 & 0 & 1 & 0 &  0 & -1\\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_s\\n_s\\k_s
\end{bmatrix}\ge\textbf{0}\\
\end{align}
Our dependence constraint system for this relation is thus (we have a dependence when $m_s=m_t$ and $k_s=n_t$):
\begin{align}
\textbf{A}_1&=\begin{bmatrix}
-1 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
-1 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & 1 & 0 & -1 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & -1 & 0 & 0 & 0\\
-1 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
-1 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
-1 & 0 & 0 & 0 & 0 & 0 & -1 & 1\\
0 & 1 & 0 & 0 & 0 & -1 & 0 & 0\\
0 & 0 & 1 & 0 & 0 & 0 & -1 & 0\\
0 & 0 & 1 & 0 & 0 & 0 &  0 & -1\\
\end{bmatrix}\\
\textbf{E}_1&=
\begin{bmatrix}
0&0&0& 1 & 0 & -1 & 0 & 0\\
0&0&0& 0 & 1 & 0 & 0 & -1\\
\end{bmatrix}\\
\textbf{0}&\le
\textbf{A}_1\begin{bmatrix}
1\\M\\N\\m_t\\n_t\\m_s\\n_s\\k_s
\end{bmatrix}\\
\textbf{0}&=\textbf{E}_1
\begin{bmatrix}
1\\M\\N\\m_t\\n_t\\m_s\\n_s\\k_s
\end{bmatrix}\\
\end{align}
However, for ease of comparison with the second example, let us shift this to use 0-based indexing:
\begin{align}
\begin{bmatrix}
0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 \\
-1 & 1 & 0 & -1 & 0 \\
-1 & 0 & 1 & 0 & -1 \\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_t\\n_t
\end{bmatrix}\ge\textbf{0}\\
\begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & -1 & 1 \\
-1 & 1 & 0 & -1 & 0 & 0 \\
-1 & 0 & 1 & 0 & -1 & 0 \\
-1 & 0 & 1 & 0 & 0 & -1 \\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_s\\n_s\\k_s
\end{bmatrix}\ge\textbf{0}\\
\end{align}
Our dependence constraint system for this relation is thus (we have a dependence when $m_s=m_t$ and $k_s=n_t$):
\begin{align}
\textbf{A}_1&=\begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
-1 & 1 & 0 & -1 & 0 & 0 & 0 & 0 \\
-1 & 0 & 1 & 0 & -1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & -1 & 1 \\
-1 & 1 & 0 & 0 & 0 & -1 & 0 & 0 \\
-1 & 0 & 1 & 0 & 0 & 0 & -1 & 0 \\
-1 & 0 & 1 & 0 & 0 & 0 & 0 & -1 \\
\end{bmatrix}\\
\textbf{E}_1&=
\begin{bmatrix}
0&0&0& 1 & 0 & -1 & 0 & 0\\
0&0&0& 0 & 1 & 0 & 0 & -1\\
\end{bmatrix}\\
\end{align}



For the second example, our individual constraint systems look like
\begin{align}
\begin{bmatrix}
0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 \\
-1 & 1 & 0 & -1 & 0 \\
-1 & 0 & 1 & 0 & -1 \\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_t\\n_t
\end{bmatrix}\ge\textbf{0}\\
\begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
-1 & 1 & 0 & -1 & 0 & 0 \\
-1 & 0 & 1 & 0 & -1 & 0 \\
-2 & 0 & 1 & 0 & -1 & -1 \\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_s\\n_s\\k_s
\end{bmatrix}\ge\textbf{0}\\
\end{align}
and our dependence constraints are
\begin{align}
\textbf{A}_2&=\begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
-1 & 1 & 0 & -1 & 0 & 0 & 0 & 0 \\
-1 & 0 & 1 & 0 & -1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
-1 & 1 & 0 & 0 & 0 & -1 & 0 & 0 \\
-1 & 0 & 1 & 0 & 0 & 0 & -1 & 0 \\
-2 & 0 & 1 & 0 & 0 & 0 & -1 & -1 \\
\end{bmatrix}\\
\textbf{E}_2&=\begin{bmatrix}
0&0&0& 1 & 0 & -1 & 0 & 0\\
-1&0&0& 0 & 1 & 0 & -1 & -1\\
\end{bmatrix}\\
\textbf{0}&\le
\textbf{A}_2
\begin{bmatrix}
1\\M\\N\\m_t\\n_t\\m_s\\n_s\\k_s
\end{bmatrix}\\
\textbf{0}&=
\textbf{E}_2\begin{bmatrix}
1\\M\\N\\m_t\\n_t\\m_s\\n_s\\k_s
\end{bmatrix}\\
\end{align}
Now, we want a schedule, i.e values of
\begin{align}
\boldsymbol{\phi}_s &= \begin{bmatrix}
\phi_{m_s}&\phi_{n_s}&\phi_{k_s}
\end{bmatrix}\\
\omega_s\\
\boldsymbol{\phi}_t &= \begin{bmatrix}
\phi_{m_t}&\phi_{n_t}
\end{bmatrix}\\
\omega_t
\end{align}
such that
\begin{align}
\phi_{m_t}m_t + \phi_{n_t}n_t + \omega_t - 
\phi_{m_s}m_s-\phi_{n_s}n_s-\phi_{k_s}k_s - \omega_s \ge 0
\end{align}
for all values of $m_t,n_t,m_s,n_s,k_s$ within the iteration space.
Farkas' lemma tell us that this is true iff we have $\lambda_0,\boldsymbol{\lambda}\ge0$ such that
\begin{align}
\lambda_0+\boldsymbol{\lambda}'\begin{bmatrix}\textbf{A}\\\textbf{E}\\-\textbf{E}\end{bmatrix}
\begin{bmatrix}1\\M\\N\\m_t\\n_t\\m_s\\n_s\\k_s\end{bmatrix}
&=
\begin{bmatrix}
\left(\omega_t-\omega_s\right)
&0&0&\boldsymbol{\phi}_t'&-\boldsymbol{\phi}_s'
\end{bmatrix}
\begin{bmatrix}1\\M\\N\\m_t\\n_t\\m_s\\n_s\\k_s\end{bmatrix}
\end{align}
Thus, if the following is satisfied
\begin{align}
\begin{bmatrix}\textbf{A}'&\textbf{E}'&-\textbf{E}'\end{bmatrix}
\boldsymbol{\lambda}
&=
\begin{bmatrix}
\left(\omega_t-\omega_s\right)-\lambda_0\\0\\0\\\boldsymbol{\phi}_t\\-\boldsymbol{\phi}_s
\end{bmatrix}
\end{align}
for some $\boldsymbol{\lambda}\ge0$, then the target iteration will take place after the source iteration for all possible values of loop induction variables (i.e. $m_t$, $k_s$, etc) and program variables (i.e. $M$ and $N$) that are within our dependence constraints.
Our focus here is on $\Delta\omega$ delta, i.e. $\omega_t-\omega_s$. The $\Delta\omega$ represents shifting across loop iterations, which is a peculiar transform we'd rather not resort to applying.
Our linear program's objective function minimizes this (as well as the $\phi$ values).

But, here we have a problem. Let's look at the first row of $\begin{bmatrix}\textbf{A}'&\textbf{E}'&-\textbf{E}'\end{bmatrix}$ for our first example:
\begin{align}
\begin{bmatrix}
1 & 0 & 0 & -1 & -1 & 0 & 0 & 0 & -1 & -1 & -1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 & 0 \\
0 & 0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & 0 & -1 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & -1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & -1 & 0 & 1 \\
\end{bmatrix}
\begin{bmatrix}
\lambda_0\\
\boldsymbol{\lambda}
\end{bmatrix}
&=
\begin{bmatrix}
\left(\omega_t-\omega_s\right)\\0\\0\\\boldsymbol{\phi}_t\\-\boldsymbol{\phi}_s
\end{bmatrix}
\end{align}
second is
\begin{align}
\begin{bmatrix}
1 & 0 & 0 & -1 & -1 & 0 & 0 & 0 & -1 & -1 & -2 & 0 & -1 & 0 & 1 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 & 0 \\
0 & 0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & 0 & -1 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & -1 & 0 & -1 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & -1 & 0 & 1 \\
\end{bmatrix}
\begin{bmatrix}
\lambda_0\\
\boldsymbol{\lambda}
\end{bmatrix}
&=
\begin{bmatrix}
\left(\omega_t-\omega_s\right)\\0\\0\\\boldsymbol{\phi}_t\\-\boldsymbol{\phi}_s
\end{bmatrix}
\end{align}

Now, one of the chief optimizations we are interested in is register tiling. As part of this, we want to hoist the loads and stores out of the inner loop. Let's say we thus want to match indices in our load and store.
So for the first example, that means selecting 
\begin{align}
\boldsymbol{\phi}_t &= \begin{bmatrix}0&1\end{bmatrix}\\
\boldsymbol{\phi}_s &= \begin{bmatrix}0&0&1\end{bmatrix}\\
\end{align}
For the second example, we have
\begin{align}
\boldsymbol{\phi}_t &= \begin{bmatrix}0&1\end{bmatrix}\\
\boldsymbol{\phi}_s &= \begin{bmatrix}0&1&1\end{bmatrix}\\
\end{align}

We could put these into a linear program solver, but it should be quickly apparent -- or at least easy to confirm -- that valid solutions are to set the $\boldsymbol{\lambda}_12=1$, and all others to $0$:
\begin{align}
\begin{bmatrix}
1 & 0 \\
0 & 0 \\
0 & 0 \\
0 & 0 \\
0 & 1 \\
0 & 0 \\
0 & 0 \\
0 & -1 \\
\end{bmatrix}\begin{bmatrix}\lambda_0\\1\end{bmatrix}
&=
\begin{bmatrix}
0\\0\\0\\0\\1\\0\\0\\-1
\end{bmatrix}\\
\begin{bmatrix}
1 & -1 \\
0 & 0 \\
0 & 0 \\
0 & 0 \\
0 & 1 \\
0 & 0 \\
0 & -1 \\
0 & -1 \\
\end{bmatrix}\begin{bmatrix}\lambda_0\\1\end{bmatrix}
&=
\begin{bmatrix}
0\\0\\0\\0\\1\\0\\-1\\-1
\end{bmatrix}\\
\end{align}
For the first example, we can set $\lambda_0=0$, while for the second we must set $\lambda_0=1$, which isn't a problem.
The $\boldsymbol{\lambda}_12$ came from the equality matrix, setting the loop induction variables into each of the array's dimensions as equal. We have two copies of these columns, with sign flipped.
This suggests matching schedules that line up with the indices used for indexing into arrays is relatively straightforward: we don't have to search for optimal $\boldsymbol{\lambda}$ values very long. If the offset is $0$, as in the first example, we can immediately find the combination that equals the corresponding schedule, without needing to run a linear program.

However, if we have an offset, the direction matters. If it is positive, $\lambda_0$ picks up the slack. Else, we need to check if some linear combination of the other rows can compensate (by running the linear program solver), possibly resorting to non-zero $\Delta\omega$.

Let's look at an example, from this very same loop nest.
Earlier, we compared the reduction load in the inner most loop with the store $A[m,n] /= U[n,n]$. This store happens after the reduction load.
But what about the store $A[m, n] = B[m, n]$? This store happens before the reduction load, thanks to being split into an earlier loop in the original source.

This means, we're in almost the same situation as before, but the source and target roles are switched; for the first example:
\begin{align}
\begin{bmatrix}
1 & 0 & 0 & -1 & -1 & 0 & 0 & 0 & -1 & -1 & -1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 & 0 \\
0 & 0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & 0 & -1 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & -1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & -1 & 0 & 1 \\
\end{bmatrix}
\begin{bmatrix}
\lambda_0\\
\boldsymbol{\lambda}
\end{bmatrix}
&=
\begin{bmatrix}
\left(\omega_t-\omega_s\right)\\0\\0\\-\boldsymbol{\phi}_s\\\boldsymbol{\phi}_t
\end{bmatrix}
\end{align}
second is
\begin{align}
\begin{bmatrix}
1 & 0 & 0 & -1 & -1 & 0 & 0 & 0 & -1 & -1 & -2 & 0 & -1 & 0 & 1 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 & 0 \\
0 & 0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & 0 & -1 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & -1 & 0 & -1 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & -1 & 0 & 1 \\
\end{bmatrix}
\begin{bmatrix}
\lambda_0\\
\boldsymbol{\lambda}
\end{bmatrix}
&=
\begin{bmatrix}
\left(\omega_t-\omega_s\right)\\0\\0\\-\boldsymbol{\phi}_s\\\boldsymbol{\phi}_t
\end{bmatrix}
\end{align}

Because we do not have an offset in the first example, the solution is exactly as it was before: $\boldsymbol{\lambda}_12=1$, and all others, including $\lambda_0$, are $0$.
However, because $\lambda$s are non-negative, our previous solution cannot be adapted.
Note that we have 
\begin{align}
\boldsymbol{\phi}_s &= \begin{bmatrix}0&1\end{bmatrix}\\
\boldsymbol{\phi}_t &= \begin{bmatrix}0&1&1\end{bmatrix}\\
\end{align}

So let's try running a linear program where we force $\Delta\omega=0$ to see if any feasible solution exists at all:
```julia
using JuMP, HiGHS
function checksat(c, C)
  # c = C * x
  model = Model(HiGHS.Optimizer)
  N = size(C,2)
  @variable(model, x[1:N] >= 0)
  @constraint(model, c .== C * x)
  @objective(model, Min, sum(x))
  optimize!(model)
  return model
end
c = [0; 0; 0; 0; -1; 0; 1; 1];
C = [
 1  0  0  -1  -1  0  0   0  -1  -1  -1   0  -1   0   1
 0  0  0   1   0  0  0   0   1   0   0   0   0   0   0
 0  0  0   0   1  0  0   0   0   1   1   0   0   0   0
 0  1  0  -1   0  0  0   0   0   0   0   1   0  -1   0
 0  0  1   0  -1  0  0   0   0   0   0   0   1   0  -1
 0  0  0   0   0  1  0   0  -1   0   0  -1   0   1   0
 0  0  0   0   0  0  1  -1   0  -1   0   0  -1   0   1
 0  0  0   0   0  0  0   1   0   0  -1   0  -1   0   1
];
checksat(c, C)
# Problem status detected on presolve: Infeasible
```
Uh oh!
So what's the problem with simply moving ahead and accepting non-zero $\Delta\omega$?
It has the effect of shifting operations across loop iterations. If we plow forward from here and satisfy all our constraints, we could end up with a solution such as:
```julia
function trisolve_ω1!(_A, _B, _U)
  A = OffsetArray(_A, -1, -1)
  B = OffsetArray(_B, -1, -1)
  U = OffsetArray(_U, -1, -1)
  M, N = size(A)
  for n = 0:N-1
    Unn = n == 0 ? zero(eltype(U)) : U[n-1, n-1]
    for m = 0:M-1
      Amn = A[m, n] = B[m, n]
      n > 0 && (A[m, n-1] /= Unn)
      for k = 0:n-1
        Amn -= A[m, k] * U[k, n]
      end
      A[m, n] = Amn
    end
  end
  let n = N
    Unn = U[n-1, n-1]
    for m = 0:M-1
      A[m, n-1] /= Unn
    end
  end
end
```
instead of the much simpler
```julia
function trisolve_ω0!(_A, _B, _U)
  A = OffsetArray(_A, -1, -1)
  B = OffsetArray(_B, -1, -1)
  U = OffsetArray(_U, -1, -1)
  M, N = size(A)
  for n = 0:N-1
    Unn = U[n, n]
    for m = 0:M-1
      Amn = B[m, n]
      for k = 0:n-1
        Amn -= A[m, k] * U[k, n]
      end
      A[m, n] = Amn / Unn
    end
  end
end
```
which requires that all $\Delta\omega=0$.

So there is obviously some failure in our analysis or representation if two representations of the same program could produce such different outcomes. It'd be great to try and understand this at a deeper level at some point. I imagine there is some insight there that I'm missing on why this should be so, that may also suggest an ideal solution.

That is, why is something as simple as redefining $k^*_t = k_t+1$, a shift of part of the polyhedra, able to solve this problem? Looking at the numbers, I can see why, but that is purely mechanical, no geometric interpretation or insight.
```julia
using LinearAlgebra
W = Matrix{Int}(I, size(C,1), size(C,1));
W[1,end] = -1;
# multiplying has the effect of subtracting the value of $k_t$
# from row 1 of `C` (the constants row)
C2 = W*C;
checksat(c, C2) # success
```

That said, we know enough to come up with solutions.
Let's first consider options.

One possibility is of course fixing it when laying out the operations. Because `trisolve_ω0!` and `trisolve_ω1!` are equivalent in order, it is of course possible to figure out how we can generate `trisolve_ω0!` from `trisolve_ω1!`'s schedule. An approach may be to check when we have non-zero omegas. If we have $\omega=1$, but no parents within the given loop iteration, we can move the operation to the end of the previous iteration. Similarly, $\omega=-1$ without children within an iteration could be shifted to the start of the next iteration.

Another possibility is to try and change the problem we solve, such that the constant offset between indices is $0$. There are a few different approaches we could imagine taking here, e.g. shifting the polyehdra, or transferring offsets into $\omega$s within the simplices we build.

The first approach is reasonable, and can be done in conjunction with the second approach.
However, I am taking the second approach first, because intuitively, we would like solving for $\Delta\omega=0$ in the LP step to also mean solving for having values associated with the same address loaded at the same time, potentially allowing us to reduce the number of load and store operations needed, or at least increase locality.

Exploring the second option in more depth, the question is, at what point do we apply the adjustment, and where?
I think we should leave the array indices, loop nest objects, and even dependence polyhedra unchanged. Instead, we can update the constraint systems to reflect applying a shift.

This is an easily invertible linear transform that can be applied to both sides of the simplices when we use them to check the direction of each dependency (we glossed over this earlier in the blog post, but the basic approach: construct mirroring simplices assuming opposite directions, and then check which is violated [we currently don't support dependencies that switch direction, which is another TODO in the future; one can break the iteration space into regions corresponding to each direction]).
We can solve this transformed problem (which can be fast - $\boldsymbol{\lambda}_12=1$, as discussed before), and then easily recover the solution in the original space
```julia-repl
julia> checksat(c, C2) # feasible, as before 

julia> csol = Int.((W.//1)\c)
8-element Vector{Int64}:
  1
  0
  0
  0
 -1
  0
  1
  1
  
julia> checksat(csol, C) # feasible
```
which of course has $\Delta\omega \ne 0$.

However, the important thing is that our dependence edges were built in the transformed problem, and it's relatively straightforward to apply these omegas as shifts in indices, without needing to do as much analysis as option $1$ implemented as a pure post-processing step.


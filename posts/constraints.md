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

However, let us consider the dependence between the store `A[m,n] = A[m,n] / U[n,n]`, and the loads used for updating in the inner loop, that is `A[m, k]` for the first version:
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

For the second example, our individual constraint systems look like
\begin{align}
\begin{bmatrix}
0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & 1\\
-1 & 1 & 0 & -1 & 0\\
-1 & 0 & 1 & 0 & -1\\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_t\\n_t
\end{bmatrix}\ge\textbf{0}\\
\begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & -1 & 1\\
-1 & 1 & 0 & -1 & 0 & 0\\
-1 & 0 & 1 & 0 & -1 & 0\\
-1 & 0 & 1 & 0 &  0 & -1\\
\end{bmatrix}\begin{bmatrix}
1\\M\\N\\m_s\\n_s\\k_s
\end{bmatrix}\ge\textbf{0}\\
\end{align}
and our dependence constraints are
\begin{align}
\textbf{A}_2&=\begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
-1 & 1 & 0 & -1 & 0 & 0 & 0 & 0\\
-1 & 0 & 1 & 0 & -1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & -1 & 1\\
-1 & 1 & 0 & 0 & 0 & -1 & 0 & 0\\
-1 & 0 & 1 & 0 & 0 & 0 & -1 & 0\\
-1 & 0 & 1 & 0 & 0 & 0 &  0 & -1\\
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
1 & -1 & -1 & 0 & 0 & -1 & -1 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
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
1 & 0 & 0 & -1 & -1 & 0 & 0 & 0 & -1 & -1 & -1 & 0 & -1 & 0 & 1 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 & 0 \\
0 & 0 & 1 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & -1 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & -1 & 0 & 0 & -1 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & -1 & 0 & -1 & 0 & 0 & -1 & 0 & 1 \\
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

Putting these into an LP solver however shows that the first is infeasible
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


```

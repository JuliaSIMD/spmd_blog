+++
title = "Farkas' Lemma"
hascode = true
date = Date(2022, 3, 24)
rss = "This post gives a worked example of Farkas Lemma applied to loop optimization"
+++

@def tags = ["farkas lemma", "polyhedral optimization", "loop-transformations"]

Consider the loop
```julia
for i = 0:I-1, j = 0:J-1
	A[i+1,j+1] = A[i+1,j] + A[i,j+1]
end
```
This loop is not trivial to paralellize, because we have dependencies on both the
`i` and `j` loop. The store into `A[i+1,j+1]` will be loaded on the next iteration
of the `j` loop by `A[i+1,j]`, and also by the next iteration of the `i` loop by
`A[i,j+1]`.

To try and optimize this loop nest, it would be helpful to have some means of modeling
the dependencies. We can use dependency polyhedra.

In general, we define polyhedra using a system of inequalities.

\begin{align}
\textbf{A}\vec{x} \le \vec{b}
\end{align}

where inequalities involving vectors apply elementwise, bold capitol letters indicate matrices,
the vector arrow $\vec{x}$ is used to indicate column vectors.

For dependence polyhedra, we define one per dependence. The aim is to describe how both
$i$ and $j$ relate in the source and target through a system of inequalities. We use the $s$ 
subscript indicate the source (parent) and $t$ subscript the target (child).
Thus, we have that the store into $A[i_s+1, j_s+1]$ is later loaded at $A[i_t+1,j_t]$,
which gives us that $i_s+1=i_t+1$ and $j_s+1=j_t$. Encoding this information as well as
the loop bounds yields the following system of inequalities, which we can interpret
as a polyhedra that gives the space of all dependent iterations (note that we can disprove
dependence by proving a depdendence polyehedra is empty, which we will demonstrate in a
future blog post). 

\begin{align}
\begin{bmatrix}
1 & 0 & 0 & 0\\
-1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & -1 & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & -1 & 0\\
0 & 0 & 0 & 1\\
0 & 0 & 0 & -1\\
1 & 0 & -1 & 0\\
-1 & 0 & 1 & 0\\
0 & 1 & 0 & -1\\
0 & -1 & 0 & 1\\
\end{bmatrix}
\begin{bmatrix}i_s\\j_s\\i_t\\j_t\end{bmatrix}
\le
\begin{bmatrix}I-1\\0\\J-1\\0\\I-1\\0\\J-1\\0\\0\\0\\-1\\1\end{bmatrix}
\end{align}

pruning a few redundant bounds, we get
\begin{align}
\begin{bmatrix}
1 & 0 & 0 & 0\\
-1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & -1 & 0 & 0\\
1 & 0 & -1 & 0\\
-1 & 0 & 1 & 0\\
0 & 1 & 0 & -1\\
0 & -1 & 0 & 1\\
\end{bmatrix}
\begin{bmatrix}i_s\\j_s\\i_t\\j_t\end{bmatrix}
\le
\begin{bmatrix}I-1\\0\\J-2\\0\\0\\0\\-1\\1\end{bmatrix}
\end{align}


Now as seen in the [orthogonalize_loops](https://spmd.org/posts/orthogonalizing_loops/) post, we like to define loop schedules as affine functions of the indices. We'll discuss more about valid 
schedules and their interpretation in a subsequent post, but for now, suffice it to say that
we need to prove things about the differences in schedules at certain levels of the loop nest.
The source iteration must occur before the target iteration, which means that we generally need
to prove one level of the loop (i.e., one row of the scheduling matrix) either equals, or
iterates before another.
Letting $\phi^l_t\left(\vec{t}\right)$ and $\phi^l_s\left(\vec{s}\right)$ be the affine schedules
for level $l$ for the target target and source respectively, such that 
\begin{align}
\phi^l_t\left(\vec{t}\right) &= 
\alpha_{t}^l + \vec{\beta}_{t}^T\vec{t}
\end{align}
and all other schedules are defined similarliy (we use the superscript $T$ is to indicate transpose).

Then we may wish to constrain ourselves to solutions that satisfy

\begin{align}
\phi^l_t\left(\vec{t}\right) - \phi^l_s\left(\vec{s}\right) \ge \delta
\end{align}
for some $\delta$.

Is there some way we can find constraints on allowed values of $\alpha_{t}$ and $\vec{\beta}_{t}$?
Enter Farkas' lemma!

According to the affine form of Farkas' lemma as given in [Uday's thesis](https://www.csa.iisc.ac.in/~udayb/publications/uday-thesis.pdf), an affine function $\psi(\vec{x})$ is non-negative everywhere in $\mathcal{D}$ iff it is a non-negative linear combination of the faces:

\begin{align}
\psi\left(\vec{x}\right)&\equiv \lambda_0 + \vec{\lambda}^T_P\left(\vec{b} - \textbf{A}\vec{x}\right),
\lambda_0\ge 0,\vec{\lambda}_P\ge 0.
\end{align}

To apply that to our problem of dependence, we let
\begin{align}
\psi\left(\left[\vec{s},\vec{t}\right]\right) 
&= \phi^l_t\left(\vec{t}\right) - \phi^l_s\left(\vec{s}\right) - \delta\\
&=
\lambda_0 + \vec{\lambda}^T_P\left(\vec{b} - \textbf{A}\vec{x}\right)
\end{align}
and then use Fourier-Motzkin elimination to remove the $\lambda$.

Returning to the first dependence of our motivating example, we have

\begin{align}
\psi^l\left(\left[\vec{s},\vec{t}\right]\right) 
&=
\psi^l\left(\left[i_s,j_s,i_t,j_t\right]\right) \\
&=
\phi^l_t\left(\left[i_t,j_t\right]\right) - \phi^l_s\left(\left[i_s,j_s\right]\right) - \delta\\
&=
\alpha_{t}^l + \beta^l_{i_t}i_t + \beta^l_{j_t}j_t
- \alpha_{s}^l - \beta_{i_s}i_s - \beta_{j_s}j_s - \delta\\
&=
\lambda_0 + 
\lambda_1\left(I-1 - i_s\right) + 
\lambda_2\left(i_s\right) + 
\lambda_3\left(J-2 - j_s\right) + 
\lambda_4\left(j_s\right) \\ &\ \ \ \ +
\lambda_5\left(i_t - i_s\right) + 
\lambda_6\left(i_s - i_t\right) + 
\lambda_7\left(-1 + j_t - j_s\right) + 
\lambda_8\left(1 + j_s - j_t\right)
\end{align}
We can produce now produce the following system of equations by matching the cofficients on
the induction variables:
\begin{align}
-\beta^l_{i_s} &= -\lambda_1 + \lambda_2 - \lambda_5 + \lambda_6\\
-\beta^l_{j_s} &= -\lambda_3 + \lambda_4 - \lambda_7 + \lambda_8\\
\beta^l_{i_t} &= \lambda_5 - \lambda_6\\
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
\alpha^l_t - \alpha^l_s - \delta &= \lambda_0 + \lambda_1\left(I-1\right) + \lambda_3\left(J-2\right)
-\lambda_7+\lambda_8
\end{align}
Along with the inequalities that $\lambda_i\ge0, i = 0,\ldots,8$.

We can perform a couple rounds of Gaussian elimination to simplify this to:
\begin{align}
\beta^l_{i_t}-\beta^l_{i_s} &= -\lambda_1 + \lambda_2\\
\beta^l_{j_t}-\beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
\beta^l_{i_t} &= \lambda_5 - \lambda_6\\
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= \lambda_0 + \lambda_1\left(I-1\right) + \lambda_3\left(J-2\right)
\end{align}
Before proceding with Fourier-Motzkin elimination to eliminate the Fourier-Motzkin multipliers (the $\lambda$s).

First to eliminate $\lambda_2$, we have the following inequalities including $\lambda_2$:
\begin{align}
\lambda_2 &\le \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
\lambda_2 &\ge \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
\lambda_2 &\ge 0
\end{align}
Eliminating $\lambda_2$ produces
\begin{align}
\beta^l_{j_t}-\beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
\beta^l_{i_t} &= \lambda_5 - \lambda_6\\
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= \lambda_0 + \lambda_1\left(I-1\right) + \lambda_3\left(J-2\right)\\
0 &\le \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
\lambda_0&\ge0\\
\lambda_1&\ge0\\
\lambda_3&\ge0\\
\lambda_4&\ge0\\
\lambda_5&\ge0\\
\lambda_6&\ge0\\
\lambda_7&\ge0\\
\lambda_8&\ge0
\end{align}
We've dropped the inequality
\begin{align}
\beta^l_{i_t}-\beta^l_{i_s} + \lambda_1 &\ge \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
\end{align}
as this simplifies to $0\ge0$. Eliminating $\lambda_4$ and $\lambda_5$, we have:
\begin{align}
\lambda_4 &\le \beta^l_{j_t}-\beta^l_{j_s} + \lambda_3\\
\lambda_4 &\ge \beta^l_{j_t}-\beta^l_{j_s} + \lambda_3\\
\lambda_4 &\ge 0\\
\lambda_5 &\le \beta^l_{i_t} + \lambda_6\\
\lambda_5 &\ge \beta^l_{i_t} + \lambda_6\\
\lambda_5 &\ge 0\\
\end{align}
yielding
\begin{align}
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= \lambda_0 + \lambda_1\left(I-1\right) + \lambda_3\left(J-2\right)\\
0 &\le \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
0 &\le \beta^l_{j_t}-\beta^l_{j_s} + \lambda_3\\
0 &\le  \beta^l_{i_t} + \lambda_6\\
\lambda_0&\ge0\\
\lambda_1&\ge0\\
\lambda_3&\ge0\\
\lambda_6&\ge0\\
\lambda_7&\ge0\\
\lambda_8&\ge0
\end{align}
Eliminating $\lambda_0$ and $\lambda_6$:
\begin{align}
\lambda_0&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_1\left(I-1\right) - \lambda_3\left(J-2\right)\\
\lambda_0&\ge\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_1\left(I-1\right) - \lambda_3\left(J-2\right)\\
\lambda_0&\ge0\\
\lambda_6 &\ge -\beta^l_{i_t}\\
\lambda_6&\ge0\\
\end{align}
Producing
\begin{align}
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
0 &\le \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
0 &\le \beta^l_{j_t}-\beta^l_{j_s} + \lambda_3\\
0 &\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_1\left(I-1\right) - \lambda_3\left(J-2\right)\\
\lambda_1&\ge0\\
\lambda_3&\ge0\\
\lambda_7&\ge0\\
\lambda_8&\ge0
\end{align}

Eliminating $\lambda_7$:
\begin{align}
\lambda_7&\le \beta^l_{j_t} + \lambda_8\\
\lambda_7&\ge \beta^l_{j_t} + \lambda_8\\
\lambda_7&\ge 0
\end{align}
Producing
\begin{align}
0 &\le \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
0 &\le \beta^l_{j_t}-\beta^l_{j_s} + \lambda_3\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_1\left(I-1\right) - \lambda_3\left(J-2\right)\\
0 &\le \beta^l_{j_t} + \lambda_8\\
\lambda_1&\ge0\\
\lambda_3&\ge0\\
\lambda_8&\ge0
\end{align}

Eliminating $\lambda_8$:
\begin{align}
\lambda_8 &\ge - \beta^l_{j_t}\\
\lambda_8&\ge0
\end{align}
it just drops out, yielding
\begin{align}
0 &\le \beta^l_{i_t}-\beta^l_{i_s} + \lambda_1\\
0 &\le \beta^l_{j_t}-\beta^l_{j_s} + \lambda_3\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_1\left(I-1\right) - \lambda_3\left(J-2\right)\\
\lambda_1&\ge0\\
\lambda_3&\ge0\\
\end{align}

Eliminating $\lambda_1$, noting that we can assume $I >= 1$, meaning either $I-1>0$ or $I-1==0$.
\begin{align}
\lambda_1\left(I-1\right) &\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_3\left(J-2\right)\\
\lambda_1 &\ge -\left(\beta^l_{i_t}+\beta^l_{i_s}\right)\\
\lambda_1&\ge0\\
\end{align}
we get
\begin{align}
0 &\le \beta^l_{j_t}-\beta^l_{j_s} + \lambda_3\\
\lambda_3&\ge0\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_3\left(J-2\right)\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - \lambda_3\left(J-2\right)
+\left(I-1\right)\left(\beta^l_{i_t}-\beta^l_{i_s}\right)
\\
\end{align}
If we assume $I-1==0$, the last equation becomes redundant. We do not need to represent this special case separately.

We can also assume that $J >= 2$, as if $J<=1$, our dependence polyhedra is empty (we have $0<=j_s<=J-2$), meaning there are no constraints that must be satisfied.
Finally, eliminating $\lambda_3$, we have
\begin{align}
\lambda_3\left(J-2\right) &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
\lambda_3\left(J-2\right) &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+\left(I-1\right)\left(\beta^l_{i_t}-\beta^l_{i_s}\right)
\\
\lambda_3&\ge0\\
\lambda_3&\ge\beta^l_{j_s}-\beta^l_{j_t}
\end{align}
yielding
\begin{align}
0&\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
0&\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+\left(I-1\right)\left(\beta^l_{i_t}-\beta^l_{i_s}\right)
\\
\left(J-2\right)\left(
\beta^l_{j_s}-\beta^l_{j_t}\right)
&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
\left(J-2\right)\left(
\beta^l_{j_s}-\beta^l_{j_t}\right)
&\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+\left(I-1\right)\left(\beta^l_{i_t}-\beta^l_{i_s}\right)
\\
\end{align}
Focusing only on this dependency, we see for example that if we want 
$\beta^l_{i_t} = \beta^l_{i_s}$ and $\beta^l_{j_t} = \beta^l_{j_s}$, then the system simplifies to
\begin{align}
\beta^l_{j_t} &\ge \delta + \alpha^l_s - \alpha^l_t
\end{align}
If this were the only dependency, e.g. if the loop was
```julia
for i = 0:I-1, j = 0:J-1
	A[i+1,j+1] = f(A[i+1,j])
end
```
this tells us that we can satisfy the dependency at level 0 via setting

\begin{align}
\alpha^0_s &= \alpha^0_t = 0\\
\beta^0_{i_s} &= \beta^0_{i_t} = 0\\
\beta^0_{j_s} &= \beta^0_{j_t} = 1\\
\end{align}
which would satisfy $\delta=1$ because $1 \ge 1 + 0 - 0$, allowing us to parallelize the inner loop (over $i$).

However, in the motivating example, we have a second dependency. We can produce its constraints similarly, referring to the second target as $u$:
\begin{align}
0&\le \beta^l_{i_u} + \alpha^l_u - \alpha^l_s - \delta\\
0&\le \beta^l_{i_u} + \alpha^l_u - \alpha^l_s - \delta
+\left(I-1\right)\left(\beta^l_{j_u}-\beta^l_{j_s}\right)
\\
\left(J-2\right)\left(
\beta^l_{i_s}-\beta^l_{i_u}\right)
&\le\beta^l_{i_u} + \alpha^l_u - \alpha^l_s - \delta\\
\left(J-2\right)\left(
\beta^l_{i_s}-\beta^l_{i_u}\right)
&\le \beta^l_{i_u} + \alpha^l_u - \alpha^l_s - \delta
+\left(I-1\right)\left(\beta^l_{j_u}-\beta^l_{j_s}\right)
\\
\end{align}
If we'd like to parallelize the inner loop via satisfying the dependency at the outer level (i.e., $\delta >= 1$ to satisfy a dependency), we can do so via the schedule
\begin{align}
\alpha^0_s &= \alpha^0_t = \alpha^0_u = 0\\
\beta^0_{i_t} &= \beta^0_{i_t} = \beta^0_{i_u} = 1\\
\beta^0_{j_s} &= \beta^0_{j_t} = \beta^0_{j_u} = 1\\
\end{align}
and then we can use
\begin{align}
\alpha^1_s &= \alpha^1_t = \alpha^1_u = 0\\
\beta^1_{i_t} &= \beta^1_{i_t} = \beta^1_{i_u} = 1\\
\beta^1_{j_s} &= \beta^1_{j_t} = \beta^1_{j_u} = 0\\
\end{align}
for the inner loop, and freely paralelize it.

Of course, note that our loop bounds change. As a final observation, note that our matrix $\boldsymbol{\beta}$ is unimodular, with determinant $-1$.
\begin{align}
\boldsymbol{\beta} &=
\begin{bmatrix}
1 & 1\\
1 & 0
\end{bmatrix}\\
\left|\boldsymbol{\beta}\right| &= -1.
\end{align}


Reference: [Effective Automatic Parallelization and Locality Optimization using the Polyehdral Model by Uday Bondhugula](https://www.csa.iisc.ac.in/~udayb/publications/uday-thesis.pdf).


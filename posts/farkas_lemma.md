+++
title = "Farkas' Lemma"
hascode = true
date = Date(2022, 3, 24)
rss = "This post gives a worked example of Farkas' Lemma applied to loop optimization"
+++

@def tags = ["farkas lemma", "polyhedral optimization", "loop-transformations"]

Consider the loop
```julia
for i = 0:I-1, j = 0:J-1
    A[i+1,j+1] = A[i+1,j] + A[i,j+1]
end
```
This loop is not trivial to parallelize because it requires serial execution of both the
`i` and `j` loops. The store into `A[i+1,j+1]` will be loaded on the next iteration
of the `j` loop by `A[i+1,j]`, and also by the next iteration of the `i` loop by
`A[i,j+1]`. We define the store into `A[i+1,j+1]` as the source, and the 
subsequent loads from this memory location as the targets which must hold the
stored value.

If the relationships between sources and targets across loop iterations
are cast as a dependency graph with coordinates corresponding to the iteration space
indices, we can model and transform those relationships using dependence polyhedra, 
allowing application of optimization techniques from linear programming and operations
research.

In general, we define polyhedra using a system of inequalities.

\begin{align}
\textbf{A}\vec{x} \le \vec{b}
\end{align}

where inequalities involving vectors apply elementwise, bold capital letters indicate matrices,
and the vector arrow $\vec{x}$ is used to indicate column vectors.

For dependence polyhedra, we define one per dependence. The aim is to describe how both
$i$ and $j$ relate in the source and target through a system of inequalities. We use the $s$ 
subscript indicate the source (parent) and $t$ subscript the target (child).
With this notation, the store into $A[i_s+1, j_s+1]$ and subsequent load from $A[i_t+1,j_t]$
imposes two constraints: $i_s+1=i_t+1$ and $j_s+1=j_t$. Encoding this information, as well as
the loop bounds, yields the following system of inequalities, which we can interpret
as a polyhedra that gives the space of all dependent iterations (note that we can disprove
dependence by proving a depdendence polyehedra is empty, which we will demonstrate in a
future blog post). 

<!-- \begin{align} -->
<!-- \begin{bmatrix} -->
<!-- 1 & 0 & 0 & 0 & -1 & 0\\ -->
<!-- -1 & 0 & 0 & 0 & 0 & 0\\ -->
<!-- 0 & 1 & 0 & 0 & 0 & -1\\ -->
<!-- 0 & -1 & 0 & 0 & 0 & 0\\ -->
<!-- 0 & 0 & 1 & 0 & -1 & 0\\ -->
<!-- 0 & 0 & -1 & 0 & 0 & 0\\ -->
<!-- 0 & 0 & 0 & 1 & 0 & -1\\ -->
<!-- 0 & 0 & 0 & -1 & 0 & 0\\ -->
<!-- 1 & 0 & -1 & 0 & 0 & 0\\ -->
<!-- -1 & 0 & 1 & 0 & 0 & 0\\ -->
<!-- 0 & 1 & 0 & -1 & 0 & 0\\ -->
<!-- 0 & -1 & 0 & 1 & 0 & 0\\ -->
<!-- \end{bmatrix} -->
<!-- \begin{bmatrix}i_s\\j_s\\i_t\\j_t\\I\\J\end{bmatrix} -->
<!-- \le -->
<!-- \begin{bmatrix}-1\\0\\-1\\0\\-1\\0\\-1\\0\\0\\0\\-1\\1\end{bmatrix} -->
<!-- \end{align} -->

<!-- pruning a few redundant bounds, we get -->
<!-- \begin{align} -->
<!-- \begin{bmatrix} -->
<!-- 1 & 0 & 0 & 0 & -1 & 0\\ -->
<!-- -1 & 0 & 0 & 0 & 0 & 0\\ -->
<!-- 0 & 1 & 0 & 0 & 0 & -1\\ -->
<!-- 0 & -1 & 0 & 0 & 0 & 0\\ -->
<!-- 1 & 0 & -1 & 0 & 0 & 0\\ -->
<!-- -1 & 0 & 1 & 0 & 0 & 0\\ -->
<!-- 0 & 1 & 0 & -1 & 0 & 0\\ -->
<!-- 0 & -1 & 0 & 1 & 0 & 0\\ -->
<!-- \end{bmatrix} -->
<!-- \begin{bmatrix}i_s\\j_s\\i_t\\j_t\\I\\J\end{bmatrix} -->
<!-- \le -->
<!-- \begin{bmatrix}-1\\0\\-2\\0\\0\\0\\-1\\1\end{bmatrix} -->
<!-- \end{align} -->

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

Pruning a few redundant bounds, we get
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


As seen in the [orthogonalize_loops](https://spmd.org/posts/orthogonalizing_loops/) post, we like to define loop schedules as affine functions of the indices. We'll discuss valid 
schedules and their interpretation in a subsequent post, but for now, suffice it to say that
we need to prove things about the differences in schedules at certain levels of the loop nest.
The source iteration must occur before the target iteration, which means that we generally need
to prove one level of the loop (i.e., one row of the scheduling matrix) either equals, or
iterates, before another.
Letting $\phi^l_t\left(\vec{t}\right)$ and $\phi^l_s\left(\vec{s}\right)$ be the affine schedules
for level $l$ for the target and source respectively, such that 
\begin{align}
\phi^l_t\left(\vec{t}\right) &= 
\alpha_{t}^l + \vec{\beta}_{t}^T\vec{t}
\end{align}
and all other schedules are defined similarly, with the superscript $T$ indicating transposition.

Then we may wish to constrain ourselves to solutions that satisfy

\begin{align}
\phi^l_t\left(\vec{t}\right) - \phi^l_s\left(\vec{s}\right) \ge \delta
\end{align}
for some $\delta$.

Is there some way we can find constraints on allowed values of $\alpha_{t}$ and $\vec{\beta}_{t}$?
Enter Farkas' lemma!

According to the affine form of Farkas' lemma given in [Uday's thesis](https://www.csa.iisc.ac.in/~udayb/publications/uday-thesis.pdf), an affine function $\psi(\vec{x})$ is non-negative everywhere in a polyhedron $\mathcal{D}$ iff it is a non-negative linear combination of the faces:

\begin{align}
\psi\left(\vec{x}\right)&\equiv \lambda_0 + \vec{\lambda}^T_P\left(\vec{b} - \textbf{A}\vec{x}\right),
\lambda_0\ge 0,\vec{\lambda}_P\ge 0.
\end{align}
That is, if some $\lambda_0\ge 0,\vec{\lambda}_P\ge 0$ exists such that the equality holds, then $\psi\left(\vec{x}\right)$ is non-negative everywhere inside the polyhedron.

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
We can now produce the following system of equations by matching 
induction variable coefficients:
\begin{align}
-\beta^l_{i_s} &= -\lambda_1 + \lambda_2 - \lambda_5 + \lambda_6\\
-\beta^l_{j_s} &= -\lambda_3 + \lambda_4 - \lambda_7 + \lambda_8\\
\beta^l_{i_t} &= \lambda_5 - \lambda_6\\
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
0 &= \lambda_1\\
0 &= \lambda_3\\
\alpha^l_t - \alpha^l_s - \delta &= \lambda_0 - \lambda_1 - 2\lambda_3 -\lambda_7+\lambda_8
\end{align}
and imposing the inequalities $\lambda_i\ge0, i = 0,\ldots,8$.

We can perform a couple rounds of Gaussian elimination to simplify this to:
\begin{align}
\beta^l_{i_t}-\beta^l_{i_s} &= \lambda_2\\
\beta^l_{j_t}-\beta^l_{j_s} &= \lambda_4\\
\beta^l_{i_t} &= \lambda_5 - \lambda_6\\
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= \lambda_0
\end{align}
before proceding with Fourier-Motzkin elimination to eliminate the Fourier-Motzkin multipliers (the $\lambda$s).

To first eliminate $\lambda_2$, we have the following inequalities including $\lambda_2$:
\begin{align}
\lambda_2 &\le \beta^l_{i_t}-\beta^l_{i_s}\\
\lambda_2 &\ge \beta^l_{i_t}-\beta^l_{i_s}\\
\lambda_2 &\ge 0
\end{align}
Eliminating $\lambda_2$ produces
\begin{align}
\beta^l_{j_t}-\beta^l_{j_s} &= \lambda_4\\
\beta^l_{i_t} &= \lambda_5 - \lambda_6\\
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= \lambda_0\\
0 &\le \beta^l_{i_t}-\beta^l_{i_s}\\
\lambda_0&\ge0\\
\lambda_4&\ge0\\
\lambda_5&\ge0\\
\lambda_6&\ge0\\
\lambda_7&\ge0\\
\lambda_8&\ge0
\end{align}
Similarly, eliminating $\lambda_0$, $\lambda_4$, $\lambda_5$, and $\lambda_7$, we have:
\begin{align}
\lambda_4 &\le \beta^l_{j_t}-\beta^l_{j_s}\\
\lambda_4 &\ge \beta^l_{j_t}-\beta^l_{j_s}\\
\lambda_4 &\ge 0\\
\lambda_5 &\le \beta^l_{i_t} + \lambda_6\\
\lambda_5 &\ge \beta^l_{i_t} + \lambda_6\\
\lambda_5 &\ge 0\\
\lambda_7&\le \beta^l_{j_t} + \lambda_8\\
\lambda_7&\ge \beta^l_{j_t} + \lambda_8\\
\lambda_7&\ge 0\\
\lambda_0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
\lambda_0 &\ge \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
\lambda_0&\ge 0
\end{align}
yielding
\begin{align}
0 &\le \beta^l_{j_t}-\beta^l_{j_s}\\
0 &\le \beta^l_{i_t} + \lambda_6\\
0 &\le \beta^l_{j_t} + \lambda_8\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
0 &\le \beta^l_{i_t}-\beta^l_{i_s}\\
\lambda_6&\ge0\\
\lambda_8&\ge0
\end{align}
Finally, eliminating $\lambda_6$ and $\lambda_8$:
\begin{align}
\lambda_6 &\ge -\beta^l_{i_t}\\
\lambda_6&\ge0\\
\lambda_8 &\ge - \beta^l_{j_t}\\
\lambda_8&\ge0
\end{align}
We realize we can just drop them.
Thus, our final equations are
\begin{align}
\beta^l_{j_s} &\le \beta^l_{j_t}\\
\beta^l_{i_s} &\le \beta^l_{i_t}\\
\delta &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s\\
\end{align}

Thus, the coefficient for the target must be at least as great as those from the source.
Additionally, if the offsets are equal, then $\beta^l_{j_t}$ must exceed our desired $\delta$.

Repeating this procedure for the other dependency would similarly yield 
\begin{align}
\beta^l_{j_s} &\le \beta^l_{j_t}\\
\beta^l_{i_s} &\le \beta^l_{i_t}\\
\delta &\le \beta^l_{i_t} + \alpha^l_t - \alpha^l_s\\
\end{align}

Note that to parallelize a hyperplane, we need $\delta$ for all unsatisfied dependencies to be 0.
One approach is to let $\beta^0_{i_t} = \beta^0_{i_s} = \beta^0_{j_t} = \beta^0_{j_s} = 1$ to satisfy both dependencies in the outer loop. Then we can freely parallelize the inner loop.

Reference: [Effective Automatic Parallelization and Locality Optimization using the Polyehdral Model by Uday Bondhugula](https://www.csa.iisc.ac.in/~udayb/publications/uday-thesis.pdf).


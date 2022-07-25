+++
title = "Ordering"
hascode = true
date = Date(2022, 7, 13)
rss = "This post describes an algorithm for answering ordering queries"
+++

@def tags = ["ordering", "linear algebra"]

Say we are given
\begin{align}
a &\le m\\
b &\le m\\
x &\ge m\\
y &\ge m
\end{align}
and want to know, is $x \ge a$? How about $x \ge y$?

A human can of course look at these equations and say "yes" and "not provable from the given information", but our goal is to develop an algorithm that can do this for arbitrary sets of inequalities, and arbitrary queries.
Our compiler will need this for some checks.

So, how can we do this?
For starters, manipulating equations tends to be easiest if we place them into matrices.
Secondarily, equalities are much easier to work with than inequalities. Thus, we introduce slack variables $s_{*}\ge0$ to turn each of these into an equality. Now, we have

\begin{align}
\begin{bmatrix}
-1 & 0 & 1 & 0 & 0 & -1 & 0 & 0 & 0\\
0 & -1 & 1 & 0 & 0 & 0 & -1 & 0 & 0\\
0 & 0 & 1 & -1 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 1 & 0 & -1 & 0 & 0 & 0 & 1\\
\end{bmatrix}
\begin{bmatrix}
a\\b\\m\\x\\y\\s_0\\s_1\\s_2\\s_3
\end{bmatrix}
&= 
\begin{bmatrix}0\\0\\0\\0
\end{bmatrix}
\end{align}

Now, say, we're given an arbitrary query with respect to our 5 variables (and possible constants). E.g., we may want to know whether $a \le x$ or $b + 1 \le y$ are known `true`.

Let us use $a \le x$, or equivalently, $x - a \ge 0$ as an example.
We could repsent this as the column vector $\begin{bmatrix}-1& 0&0&1&0\end{bmatrix}^T$.

So, to be able to answer queries with respect to $x$ and $a$, the first thing we'd need to be able to do isolate them in an expression.
That means, we'll need to take some linear combination of the constraints to get the vector $\textbf{q}=\begin{bmatrix}-1& 0&0&1&0\end{bmatrix}^T$.

But first, we should put our constraints into Hermite Normal Form and then drop all the zeroed constraints, or employ some other similar technique, in order to make sure that all remaining constraints are linearly independent. Guaranteeing they're linearly independent may simplify work later on.
In the example above, we already have full row rank (4), so we'll skip this step in this blog post.

Now, to isolate our query, we have the problem

\begin{align}
\begin{bmatrix}
-1 & 0 & 0 & 0\\
0 & -1 & 0 & 0\\
1 & 1 & 1 & 1\\
0 & 0 & -1 & 0\\
0 & 0 & 0 & -1\\
-1 & 0 & 0 & 0\\
0 & -1 & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1\\
\end{bmatrix}
\textbf{x}&=
\begin{bmatrix}\textbf{q}\\\textbf{z}\end{bmatrix}
\end{align}
So we're trying to find a linear combination $\textbf{x}$ of constraints to isolate our query $\textbf{q}$, but then what shall $\textbf{z}$ be?

$\textbf{z}$ corresponds to our four slack variables.

Now, to answer queries of the $\ge$ variety, we want all slack variables to be $\le0$.
That is if we have $w - s_{*} = d$, then $w\ge d$.

We can focus on only $\ge$, as the other queries of interest can be reframed as $\ge$, e.g. by adjusting constant offsets or flipping signs.

So, our current idea for proving (or failing to prove) a query is to find the linear combination of existing constraints that produces our query, and then check if all slack variables are negative.

To assist us in this quest, we add augment colums and set $\textbf{z}=\textbf{0}$:
\begin{align}
\begin{bmatrix}
-1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & -1 & 0 & 0 & 0 & 0 & 0 & 0\\
1 & 1 & 1 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & -1 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & -1 & 0 & 0 & 0 & 0\\
-1 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & -1 & 0 & 0 & 0 & 1 & 0 & 0\\
0 & 0 & 1 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 1\\
\end{bmatrix}
\textbf{x}&=
\begin{bmatrix}\textbf{q}\\0\\0\\0\\0\end{bmatrix}
\end{align}

Now, the linear combination of columns $\textbf{x}$ must produce a value of $0$ for each slack variable. We can divide this into both the part contributing to our query, and our augment. Letting $V$ be the total number of variables, and $E$ the total number of equations, using 0-based indexing, and referring to the above expression more succinctly as $\textbf{Ax}=\textbf{q}$:

\begin{align}
0 &= s_{e}\\
&= \sum_{c=0}^{E-1} A_{e+V,c} x_c +
\sum_{c=0}^{E-1} A_{e+V,c+E} x_c\\
&= x_e + \sum_{c=0}^{E-1} A_{e+V,c} x_c\\
x_e &= -\sum_{c=0}^{E-1} A_{e+V,c} x_c\\
\end{align}

Note that $\sum_{c=0}^{E-1} A_{e+V,c} x_c$ is precisely the quantity we require to be $\le 0$ for our query to be true; it is the value of the slack variable in the linear combination of constraints that produces our query.
So, to check if this is satisfied, all we must do is check that $x_e\ge0$.

Okay, but how do we solve for $x$? If $\textbf{A}$ has full column rank, this is trivial: we either have a solution, or no solutions. If we have no solutions, we obviously fail. If we have one solution, it is trivial to check whether all $x_e\ge0, e = E,\ldots,2E-1$.

But what if $\textbf{A}$ does not have full column rank?
Again, we have the problem $\textbf{Ax}=\textbf{q}$. Our first step will be to ensure that it has full row rank, which we can do by multiplying by some matrix $\textbf{U}$ if necessary, e.g., calculated from Hermite Normal Form, discarding zeroed rows. Letting $\textbf{UA}=\textbf{H}$, and $\textbf{Uq}=\textbf{b}$, we now have $\textbf{Hx}=\textbf{b}$.

Letting $M$ and $N$ be the number of rows and columns of $\textbf{H}$, respectively. Note that $M$ is also the rank of $\textbf{H}$, as we have guaranteed that it has full row rank. Because the column rank is deficient, we know $N > M$.

We can construct a non-singular (but not necessarilly unimodulary) $N\times N$ integer matrix $\textbf{V}$, where the first $M$ columns diagonalize $\textbf{H}$, while the remaining $N-M$ columns are the nullspace of $\textbf{H}$. Calling these two blocks $\textbf{V}_1$ and $\textbf{V}_2$, we have $\textbf{V} = \begin{bmatrix}\textbf{V}_1&\textbf{V}_2\end{bmatrix}$.

With this, we have

\begin{align}
\textbf{b} &=
\textbf{HVV}^{-1}\textbf{x}\\
&=
\textbf{H}
\begin{bmatrix}\textbf{V}_1&\textbf{V}_2\end{bmatrix}
\textbf{V}^{-1}\textbf{x}\\
&=
\begin{bmatrix}\textbf{D}&\textbf{0}_{N\times \left(N-M\right)}\end{bmatrix}
\textbf{V}^{-1}\textbf{x}\\
\end{align}
Where $\textbf{D}$ is diagonal.

Now, let

\begin{align}
\textbf{y}_1 &= \textbf{D}^{-1}\textbf{b}\\
\textbf{V}^{-1}\textbf{x} &= \begin{bmatrix}\textbf{y}_1\\\textbf{y}_2\end{bmatrix}\\
\textbf{b} &=
\begin{bmatrix}\textbf{D}&\textbf{0}_{N\times \left(N-M\right)}\end{bmatrix}
\begin{bmatrix}\textbf{D}^{-1}\textbf{b}\\\textbf{y}_2\end{bmatrix}\\
&= \textbf{b} + \textbf{0}\\
&= \textbf{b}
\end{align}

Therefore, letting 

\begin{equation}
\textbf{V}^{-1}\textbf{x} = \begin{bmatrix}\textbf{D}^{-1}\textbf{b}\\\textbf{y}_2\end{bmatrix}\\
\end{equation}

is a valid solution for all values of $\textbf{y}_2$. Now, solving for $x$, as it is $x$ that we need to prove our queries, we have

\begin{align}
\textbf{V}^{-1}\textbf{x} &= \begin{bmatrix}\textbf{D}^{-1}\textbf{b}\\\textbf{y}_2\end{bmatrix}\\
\textbf{x}
&=
\begin{bmatrix}\textbf{V}_1&\textbf{V}_2\end{bmatrix}
\begin{bmatrix}\textbf{D}^{-1}\textbf{b}\\\textbf{y}_2\end{bmatrix}\\
\textbf{x}
&=
\textbf{V}_1\textbf{D}^{-1}\textbf{b} +
\textbf{V}_2\textbf{y}_2\\
\end{align}
As $\textbf{V}_2$ is the nullspace of $\textbf{H}$, $\textbf{y}_2$ can moves through the nullspace of $\textbf{H}$.

Now, recognizing that we need the last $E$ entries of $\textbf{x}$ to be $\ge 0$ but don't care about the first $V$ elements, we left-multiply by

\begin{equation}
\textbf{J} = \begin{bmatrix}\textbf{0}_{E\times E}&\textbf{I}_{E\times E}\end{bmatrix}
\end{equation}

to discard the elements of $\textbf{x}$ we aren't interested in. For simplicity, let

\begin{align}
\textbf{x}^* &= \textbf{Jx}\\
\textbf{c} &= \textbf{JV}_1\textbf{D}^{-1}\textbf{b}\\
\textbf{W} &= -\textbf{JV}_2\\
\textbf{x}^* &= c - \textbf{Wy}_2
\end{align}

Now, as what we need is $\textbf{x}^*\ge 0$

\begin{align}
\textbf{0} &\le c - \textbf{Wy}_2\\
\textbf{Wy}_2 &\le c\\
\end{align}


So, this reduces the problem to one of simply checking if $\textbf{Wy}_2 \le c$ is satisfiable, where $\textbf{W}$ is known from the initial constraints, and $\textbf{c}$ is computed from the query (and most of the computation can be done once on initialization; for a particular query, all we need is a matrix-vector multiply). Checking feasibility is the first step to solving any linear program, but we can write specialized software for polyhedra that handles the unbounded nature of $\textbf{y}_2$ without needing to duplicate the elements of $\textbf{W}$ in a tableau.



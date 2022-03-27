+++
title = "Farkas' Coninued"
hascode = true
date = Date(2022, 3, 24)
rss = "Discussion of Farkas Lemma with symbolic variables on either side-of-zero"
+++

@def tags = ["farkas lemma", "polyhedral optimization", "loop-transformations"]


This is a slightly modification of the example from the [previous Farkas post](https://spmd.org/posts/farkas_lemma/), adjusting the loop bounds to something with less obvious sign:
```julia
for i = Lᵢ:Uᵢ, j = Lⱼ:Uⱼ
    A[i+1,j+1] = A[i+1,j] + A[i,j+1]
end
```
This gives us the following dependence polyhedra for the first of the dependences, going from 
the store into `A[i+1,j+1]` into the load `A[i+1,j]`:
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
\begin{bmatrix}U_i\\-L_i\\U_j-1\\-L_j\\0\\0\\-1\\1\end{bmatrix}
\end{align}

In the previous example, the lower bounds were `0`, and the upper bounds were thus necessarilly non-negative -- if they were negative, our loop wouldn't iterate, and our dependence polyhedra would've been empty, meaning whatever transforms we applied would be irrelevant.

Here, however, it is only required that $L_i \le U_i$ and $L_j \le U_j$.

Why do we care about the sign? When performing Fourier-Motzkin, we need to know the sign of a coefficient to know whether the given bound is a lower or an upper bound!

However, intuitively, it seems like this should not really change anything from the previous example. Hence, why I am working through an example to try and build some intuition.

We start with applying Farkas lemma, yielind:
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
\lambda_1\left(U_i - i_s\right) + 
\lambda_2\left(i_s - L_i\right) + 
\lambda_3\left(U_j - j_s\right) + 
\lambda_4\left(j_s - L_j\right) \\ &\ \ \ \ +
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
\alpha^l_t - \alpha^l_s - \delta &= 
\lambda_0 + U_i\lambda_1 - L_i\lambda_2 + U_j\lambda_3 - L_j\lambda_4
-\lambda_7+\lambda_8
\end{align}
We can start off as we did before, first applying Gaussian elmination:
\begin{align}
\beta^l_{i_t} - \beta^l_{i_s} &= -\lambda_1 + \lambda_2\\
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
\beta^l_{i_t} &= \lambda_5 - \lambda_6\\
\beta^l_{j_t} &= \lambda_7 - \lambda_8\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= 
\lambda_0 + U_i\lambda_1 - L_i\lambda_2 + U_j\lambda_3 - L_j\lambda_4
\end{align}
As each of $\lambda_1$, $\lambda_2$, $\lambda_3$, and $\lambda_4$ have coefficients with unknown signs, we will first remove the other 5 farkas multipliers to simplify the problem.

Removing $\lambda_5$ and $\lambda_7$:
\begin{align}
\lambda_5 &\le \beta^l_{i_t} + \lambda_6\\
\lambda_5 &\ge \beta^l_{i_t} + \lambda_6\\
\lambda_5 &\ge 0\\
\lambda_7 &\le \beta^l_{j_t} + \lambda_8\\
\lambda_7 &\ge \beta^l_{j_t} + \lambda_8\\
\lambda_7 &\ge 0
\end{align}
We have:
\begin{align}
\beta^l_{i_t} - \beta^l_{i_s} &= -\lambda_1 + \lambda_2\\
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= 
\lambda_0 + U_i\lambda_1 - L_i\lambda_2 + U_j\lambda_3 - L_j\lambda_4\\
0 &\le \beta^l_{i_t} + \lambda_6\\
0 &\le \beta^l_{j_t} + \lambda_8\\
\lambda_0 &\ge 0\\
\lambda_1 &\ge 0\\
\lambda_2 &\ge 0\\
\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\lambda_6 &\ge 0\\
\lambda_8 &\ge 0\\
\end{align}
The inequalities with $\lambda_6$ and $\lambda_8$:
\begin{align}
\lambda_6 &\ge -\beta^l_{i_t}\\
\lambda_8 &\ge -\beta^l_{j_t}\\
\lambda_6 &\ge 0\\
\lambda_8 &\ge 0\\
\end{align}
So we can simply drop these, as before, producing:
\begin{align}
\beta^l_{i_t} - \beta^l_{i_s} &= -\lambda_1 + \lambda_2\\
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta &= 
\lambda_0 + U_i\lambda_1 - L_i\lambda_2 + U_j\lambda_3 - L_j\lambda_4\\
\lambda_0 &\ge 0\\
\lambda_1 &\ge 0\\
\lambda_2 &\ge 0\\
\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
Now, the inequalities with $\lambda_0$:
\begin{align}
\lambda_0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 + L_i\lambda_2 - U_j\lambda_3 + L_j\lambda_4\\
\lambda_0 &\ge \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 + L_i\lambda_2 - U_j\lambda_3 + L_j\lambda_4\\
\lambda_0 &\ge 0\\
\end{align}
Leaving us with only the $\lambda$s featuring coefficients of unknown sign:

\begin{align}
\beta^l_{i_t} - \beta^l_{i_s} &= -\lambda_1 + \lambda_2\\
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 + L_i\lambda_2 - U_j\lambda_3 + L_j\lambda_4\\
\lambda_1 &\ge 0\\
\lambda_2 &\ge 0\\
\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}


The obvious approach here is to try all $3^4$ possiblities, where $4$ is the number of $\lambda$s, and $3$ for each of negative, 0, and positive. Of course, we can prune this almost immediately, by noting that if a lower bound is 0 or positive, this restricts the possibilities of the upper bound's sign.

Let's take a "breadth-first" approach to see how things evolve, elimating $\lambda_2$ under each of the assumptions. Solving for $\lambda_2$, we get
\begin{align}
-L_i\lambda_2 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 - U_j\lambda_3 + L_j\lambda_4\\
\lambda_2 &\le \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\\
\lambda_2 &\ge \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\\
\lambda_2 &\ge 0\\
\end{align}
Assuming $L_i < 0$, we get
\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 - U_j\lambda_3 + L_j\lambda_4\\
0 &\le \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\\
-L_i\left(\beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\right)
&\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 - U_j\lambda_3 + L_j\lambda_4\\
\lambda_1 &\ge 0\\
\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
Assuming $L_i = 0$:
\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 - U_j\lambda_3 + L_j\lambda_4\\
0 &\le \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\\
\lambda_1 &\ge 0\\
\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
Assuming $L_i > 0$:
\begin{align}
\lambda_2 &\le \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\\
L_i\lambda_2 &\ge -\beta^l_{j_t} - \alpha^l_t + \alpha^l_s + \delta + U_i\lambda_1 + U_j\lambda_3 - L_j\lambda_4\\
\lambda_2 &\ge \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\\
\lambda_2 &\ge 0\\
\end{align}

\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
0&\le \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\\
-\beta^l_{j_t} - \alpha^l_t + \alpha^l_s + \delta + U_i\lambda_1 + U_j\lambda_3 - L_j\lambda_4
&\le L_i\left( \beta^l_{i_t} - \beta^l_{i_s} + \lambda_1\right)\\
\lambda_1 &\ge 0\\
\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
All three branches here are nearly the same equation, but we have an additional
\begin{align}
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_i\lambda_1 - U_j\lambda_3 + L_j\lambda_4\\
\end{align}
in the $L_i < 0$ case.
To try and tap down on the potential exponential growth of considering all possibilities, we'll try taking the intersection of all inequalities. In this case, that is equivalent with using those from the $L_i<0$ case.

We can observe that if we're matching a coefficient with unknown sign vs an equality, we'll get the same outcome either way.

Additionally, solving for $\lambda_1$, we now get
\begin{align}
\left(U_i-L_i\right)\lambda_1&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
U_i\lambda_1 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4\\
\lambda_1 &\ge \beta^l_{i_s} - \beta^l_{i_t}\\
\lambda_1 &\ge 0\\
\end{align}
Looking at these equations, we can also see that it makes sense: if $L_i < 0$, then 
\begin{align}
\left(U_i-L_i\right)\lambda_1 > U_i\lambda_1
\end{align}
while if $L_i > 0$, we have
\begin{align}
\left(U_i-L_i\right)\lambda_1 < U_i\lambda_1
\end{align}
so, depending on what we assume about $L_i$, we'd actually be able to prune one of these two bounds.
With the lack of information, we keep both.

Now, we also do not know $U_i$'s value with respect to $0$, except that $U_i \ge L_i$ (meaning if $L_i > 0$, then $U_i > 0$).
So, now, let us again split across these possibilities; if $U_i < 0$:
\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
\left(U_i-L_i\right)\left(\beta^l_{i_s} - \beta^l_{i_t}\right)
&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

\left(L_i-U_i\right)
\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4\right]
&\le -U_i
\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\right]\\

\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
Simplifying...
\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0
&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
U_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0&\le
-L_i
\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4
+ U_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right]
\\


\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
Things are starting to look a little hairy. Neither the second nor third bounds above render the other redundant, even though we know $U_i > L_i$, because we don't know the sign of $\left( \beta^l_{i_t} - \beta^l_{i_s} \right)$.
The fourth row, on the other hand, is redundant. Because $L_i < U_i < 0$, we have $-L_i>0$. Thus, we can divide both sides of the inequality by $-L_i$, dropping the fourth as equivalent to the third.
Meaning, we have
\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\
0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0
&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
U_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
as the final set of bounds under the $U_i < 0$ condition.

If $U_i > 0$, we'd instead have
\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\

0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4\\
\left(\beta^l_{i_s} - \beta^l_{i_t}\right)U_i &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4\\
\left(\beta^l_{i_s} - \beta^l_{i_t}\right)\left(U_i-L_i\right) &\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}
which we can simplify to (note that the fourth and fifth bounds above are equivalent):
\begin{align}
\beta^l_{j_t} - \beta^l_{j_s} &= -\lambda_3 + \lambda_4\\

0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4\\
0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4 + U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\

\lambda_3 &\ge 0\\
\lambda_4 &\ge 0\\
\end{align}

This is the same as in the $U_i<0$ case, but with the addition of
\begin{align}
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta - U_j\lambda_3 + L_j\lambda_4\\
\end{align}
The $U_i=0$ case does not add anything.


Pressing on, we take the intersection (equivalent to taking the bounds we got from $U_i>0$), solving for $\lambda_3$, mostly to mix things up and pick the Farkas multiplier associated with an upper bound.
\begin{align}
\lambda_3  &\le \lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)\\
\lambda_3  &\ge \lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)\\

U_j\lambda_3 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4\\
U_j\lambda_3 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
U_j\lambda_3 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 + U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
\lambda_3 &\ge 0\\
\end{align}

Assuming $U_j < 0$...

\begin{align}
0&\le \lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)\\

-\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4
\right] &\le
-U_j\left[\lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\right]\\
-\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\right] &\le
-U_j\left[\lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\right]\\
-\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 + U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)
\right] &\le
-U_j\left[\lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\right]\\
\lambda_4 &\ge 0\\
\end{align}
Simplifying...
\begin{align}

\left(U_j-L_j\right)\lambda_4
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)\\
\left(U_j-L_j\right)\lambda_4
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
\left(U_j-L_j\right)\lambda_4
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
\lambda_4 &\ge\beta^l_{j_t} - \beta^l_{j_s}\\
\lambda_4 &\ge 0\\
\end{align}

Assuming $U_j > 0$...
\begin{align}
U_j\left[\lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)\right]&\le
 \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4\\
U_j\left[\lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)\right]&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
U_j\left[\lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)\right]&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 + U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
0 &\le\lambda_4 - \left(\beta^l_{j_t} - \beta^l_{j_s}\right)\\
0 &\le \beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4\\
0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0 &\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + L_j\lambda_4 + U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
\lambda_4 &\ge 0\\
\end{align}
Simplifying...
\begin{align}
\left(U_j-L_j\right)\lambda_4
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\

\left(U_j-L_j\right)\lambda_4
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

\left(U_j-L_j\right)\lambda_4
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
\lambda_4&\ge\beta^l_{j_t} - \beta^l_{j_s}\\
\lambda_4 &\ge 0\\

L_j\lambda_4&\ge
-\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\right]
\\
L_j\lambda_4&\ge
-\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right]
\\
L_j\lambda_4&\ge
-\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\right]
\\
\end{align}

We can again take the intersection, which gives us the same equations as assuming $U_j>0$.

So, let us now assume $L_j<0$:
\begin{align}
\left(U_j-L_j\right)\left(\beta^l_{j_t} - \beta^l_{j_s}\right)&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)\\

\left(U_j-L_j\right)\left(\beta^l_{j_t} - \beta^l_{j_s}\right)&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

\left(U_j-L_j\right)\left(\beta^l_{j_t} - \beta^l_{j_s}\right)&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\

-L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\

-L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

-L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\

0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
\end{align}
Simplifying (we're able to drop 3 redundant bounds)...
\begin{align}
0&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)\\

0&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

0&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\

0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\
0&\le\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\
\end{align}

So, let us now assume $L_j>0$:
\begin{align}
0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\

0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\

\left(U_j-L_j\right)\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\

\left(U_j-L_j\right)\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

\left(U_j-L_j\right)\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\


\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\right]\\

\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\right]\\
\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)
\right]\\


\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\right]\\

\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\right]\\

\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)
\right]\\


\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\right]\\

\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\right]\\

\left(L_j-U_j\right)\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\right]
&\le
L_j\left[
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)
\right]\\
\end{align}
Simplifying...
\begin{align}
0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\

0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\

0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta
+
L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\

0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\\

0
&\le
\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
L_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\\

0
&\le
U_j\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\right]
+
L_jU_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\
0
&\le
U_j\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\right]
+
L_jU_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
+
L_jL_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\\
0
&\le
U_j\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta\right]
+
L_jU_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
+
L_jU_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)
\\

0
&\le
U_j\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right]
+
L_jU_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
-
L_jL_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\\

0
&\le
U_j\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right] +
L_jU_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\

0
&\le
U_j
\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)\right]
+
L_j\left[
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right) -
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right)
\right]\\


0
&\le
U_j\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\right]
+
L_j\left[
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)-
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)
\right]\\

0
&\le
U_j
\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\right]
+
L_j\left[
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right) +
L_i\left( \beta^l_{i_t} - \beta^l_{i_s} \right) -
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)
\right]\\

0
&\le
U_j
\left[\beta^l_{j_t} + \alpha^l_t - \alpha^l_s - \delta + 
U_i\left(\beta^l_{i_t} - \beta^l_{i_s}\right)\right]
+
L_j
U_j\left(\beta^l_{j_t} - \beta^l_{j_s}\right)
\\
\end{align}
This is a bit of a mess.
Hopefully I messed up my arithmetic somewhere, otherwise this is really starting to look exponential.
If we restrict ourselves to $\beta^l_{j_t} = \beta^l_{j_s}$ and $\beta^l_{i_t} = \beta^l_{i_s}$, much of this drops out, but we do not want to restrict ourselves to this in general.

This demonstrates the importance of canonicalizing loops in the internal representation to set $L_i=L_j=0$.



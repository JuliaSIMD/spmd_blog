+++
title = "Unimodular Matrices and Loops"
hascode = true
date = Date(2022, 3, 2)
rss = "The significance of unimodular matrices in the context of loop optimizations"
+++

@def tags = ["unimodular", "loop-transformation"]

A unimodular matrix is a square integer matrix with determinant $\pm 1$. With
some basic linear algebra lemmas, we can show that a square integer matrix is
unimodular if and only if its inverse exists and the inverse matrix is also an
integer matrix.

Let $A$ be an unimodular matrix. As the determinant of an integer matrix must
also be an integer, $\adj(A)$ is also an integer matrix. It follows that $A^{-1}
= \frac{1}{\det{A}} \adj(A)$ is an integer matrix.

To show the converse, let $A$ be an integer square matrix such that $A^{-1}$ is
an integer matrix as well. Since $\det(A) \det(A^{-1}) = \det(A A^{-1}) = \det(I)
= 1$ and $\det(A), \det(A^{-1}) \in \ZZ$, we can conclude that $A$ and $A^{-1}$
are unimodular.


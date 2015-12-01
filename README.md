# LinearAlgebra
Basic dense linear algebra routines and classes optimized for ease of use
- Matrices and Vectors
- Abstracted Matrix class to allow for possible Sparse Matrix extension
- Get a basis for a vector space via QR (Gram-Schmidt, Householder, Givens, Randomized method)
- Compute LU, SVD, or Cholesky factorization
- Perform eigenvalue decomposition
- Solve linear systems:
  - Direct Methods: Forward/Back substitution, Diagonal solve, Tridiagonal solve, Unitary solve
  - Iterative Methods: Steepest Descent, Conjugate Gradient, GMRES, Jacobi, Lanczos, BCG, QMR, Arnoldi, Bi-CGSTAB, TFQMR

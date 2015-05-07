/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ASymGShf.cc.
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in buckling mode using the AREig
   function.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in buckling
      mode, where A and B are obtained from the finite element
      discretization of the 1-dimensional discrete Laplacian
                                  d^2u / dx^2
      on the interval [0,1] with zero Dirichlet boundary conditions
      using piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: upper triangular part of matrix A
                                  stored in CSC format.
      {nnzB, irowB, pcolB, valB}: upper triangular part of matrix B
                                  stored in CSC format.

   3) Library called by this example:

      The SuperLU package is called by AREig to solve some linear
      systems involving (A-sigma*B).

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lsmatrxc.h       SymmetricMatrixC, a function that generates
                       matrix A in CSC format.
      lsmatrxd.h       SymmetricMatrixD, a function that generates
                       matrix B in CSC format.
      areig.h          The AREig function definition.
      asymsol.h        The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "lsmatrxc.h"
#include "lsmatrxd.h"
#include "areig.h"
#include "asymsol.h"


int main()
{

  // Defining variables;

  int    n;              // Dimension of the problem.
  int    nconv;          // Number of "converged" eigenvalues.
  int    nnzA,   nnzB;   // Number of nonzero elements in A and B.
  int    *irowA, *irowB; // pointer to an array that stores the row
                         // indices of the nonzeros in A and B.
  int    *pcolA, *pcolB; // pointer to an array of pointers to the
                         // beginning of each column of A (B) in valA (valB).
  double *valA,  *valB;  // pointer to an array that stores the nonzero
                         // elements of A and B.
  double EigVal[101];    // Eigenvalues.
  double EigVec[1001];   // Eigenvectors stored sequentially.
  char   uplo;           // Variable that indicates whether the upper
                         // (uplo='U') ot the lower (uplo='L') part of
                         // A and B will be supplied to AREig.

  // Creating matrices A and B.

  n = 100;
  uplo = 'U';
  SymmetricMatrixC(n, nnzA, valA, irowA, pcolA, uplo);
  SymmetricMatrixD(n, nnzB, valB, irowB, pcolB, uplo);

  // Finding the four eigenvalues of A nearest to 1.0 and the
  // related eigenvectors.

  nconv = AREig(EigVal, EigVec, n, nnzA, valA, irowA, pcolA,
                nnzB, valB, irowB, pcolB, uplo, 'B', 1.0, 4);

  // Printing solution.

  Solution(nconv, n, nnzA, valA, irowA, pcolA, nnzB,
           valB, irowB, pcolB, uplo, EigVal, EigVec);

} // main.


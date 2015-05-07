/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ASymReg.cc.
   Example program that illustrates how to solve a real
   symmetric standard eigenvalue problem in regular mode
   using the AREig function.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular 
      mode, where A is derived from the standard central difference
      discretization of the 2-dimensional Laplacian on the unit 
      square with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      {nnz, irow, pcol, A}: upper triangular part of matrix A 
                            stored in CSC format.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lsmatrxa.h       SymmetricMatrixA, a function that generates
                       matrix A in CSC format.
      areig.h          The AREig function definition.
      asymsol.h        The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "lsmatrxa.h"
#include "areig.h"
#include "asymsol.h"


int main()
{

  // Defining variables;

  int     nx;
  int     n;           // Dimension of the problem.
  int     nconv;       // Number of "converged" eigenvalues.
  int     nnz;         // Number of nonzero elements in A.
  int*    irow;        // pointer to an array that stores the row
                       // indices of the nonzeros in A.
  int*    pcol;        // pointer to an array of pointers to the
                       // beginning of each column of A in vector A.
  double* A;           // pointer to an array that stores the
                       // nonzero elements of A.
  double EigVal[101];  // Eigenvalues.
  double EigVec[1001]; // Eigenvectors stored sequentially.
  char    uplo;        // Variable that indicates whether the upper
                       // (uplo='U') ot the lower (uplo='L') part of
                       // A will be stored in A, irow and pcol.

  // Creating a 100x100 matrix.

  nx = 10;
  uplo = 'U';
  SymmetricMatrixA(nx, n, nnz, A, irow, pcol, uplo);

  // Finding the four eigenvalues with smallest magnitude and
  // the related eigenvectors.

  nconv = AREig(EigVal, EigVec, n, nnz, A, irow, pcol, uplo, 4, "SM");

  // Printing solution.

  Solution(nconv, n, nnz, A, irow, pcol, uplo, EigVal, EigVec);

} // main.


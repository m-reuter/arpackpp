/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ANSymReg.cc.
   Example program that illustrates how to solve a real
   nonsymmetric standard eigenvalue problem in regular mode
   using the AREig function.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is derived from the standard central difference
      discretization of the 2-dimensional convection-diffusion operator
                       (Laplacian u) + rho*(du/dx)
      on a unit square with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lnmatrxb.h       BlockTridMatrix, a function that generates
                       matrix A in CSC format.
      areig.h          The AREig function definition.
      ansymsol.h       The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "lnmatrxb.h"
#include "areig.h"
#include "ansymsol.h"


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
  double EigValR[101]; // Real part of the eigenvalues.
  double EigValI[101]; // Imaginary part of the eigenvalues.
  double EigVec[1001]; // Eigenvectors stored sequentially.

  // Creating a double precision 100x100 matrix.

  nx = 10;
  BlockTridMatrix(nx, n, nnz, A, irow, pcol);

  // Finding the four eigenvalues with largest magnitude and 
  // the related eigenvectors.

  nconv = AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol, 4);

  // Printing solution.

  Solution(nconv, n, nnz, A, irow, pcol, EigValR, EigValI, EigVec);

} // main.

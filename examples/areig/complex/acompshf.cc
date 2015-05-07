/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ACompShf.cc
   Example program that illustrates how to solve a complex standard
   eigenvalue problem in shift and invert mode using the AREig function.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and invert
      mode, where A is derived from the central difference discretization
      of the 1-dimensional convection-diffusion operator
                        (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by AREig to solve some linear 
      systems involving (A-sigma*I). This is needed to implement 
      the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      lcmatrxb.h       CompMatricB, a function that generates matrix
                       A in CSC format.
      areig.h          The AREig function definition.
      acompsol.h       The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arcomp.h"
#include "lcmatrxb.h"
#include "areig.h"
#include "acompsol.h"

int main()
{

  // Defining variables;

  int               n;            // Dimension of the problem.
  int               nnz;          // Number of nonzero elements in A.
  int               nconv;        // Number of "converged" eigenvalues.
  int*              irow;         // pointer to an array that stores the row
                                  // indices of the nonzeros in A.
  int*              pcol;         // pointer to an array of pointers to the
                                  // beginning of each column in vector A.
  arcomplex<double> rho;          // Parameter used by CompMatrixB.
  arcomplex<double> *A;           // pointer to an array that stores the
                                  // nonzero elements of A.
  arcomplex<double> EigVal[101];  // Eigenvalues.
  arcomplex<double> EigVec[1001]; // Eigenvectors stored sequentially.

  // Creating a complex matrix.

  n   = 100;
  rho = 10.0;
  CompMatrixB(n, rho, nnz, A, irow, pcol);

  // Finding the four eigenvalues of A nearest to 0.0 and the 
  // related eigenvectors.

  nconv = AREig(EigVal, EigVec, n, nnz, A, irow, pcol,
                arcomplex<double>(0.0, 0.0), 4L);

  // Printing solution.

  Solution(nconv, n, nnz, A, irow, pcol, EigVal, EigVec);

} // main.

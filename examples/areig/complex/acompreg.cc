/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ACompReg.cc.
   Example program that illustrates how to solve a complex standard
   eigenvalue problem in regular mode using the AREig function.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is obtained from the standard central difference
      discretization of the convection-diffusion operator
                    (Laplacian u) + rho*(du / dx)
      on the unit squre [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.

   3) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      lcmatrxa.h       CompMatrixA, a function that generates matrix
                       A in CSC format.
      areig.h          The AREig function definition.
      acompsol.h       The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arcomp.h"
#include "areig.h"
#include "lcmatrxa.h"
#include "acompsol.h"


int main()
{

  // Defining variables;

  int               nx;
  int               n;            // Dimension of the problem.
  int               nnz;          // Number of nonzero elements in A.
  int               nconv;        // Number of "converged" eigenvalues.
  int*              irow;         // pointer to an array that stores the row
                                  // indices of the nonzeros in A.
  int*              pcol;         // pointer to an array of pointers to the
                                  // beginning of each column in vector A.
  arcomplex<double> *A;           // pointer to an array that stores the
                                  // nonzero elements of A.
  arcomplex<double> EigVal[101];  // Eigenvalues.
  arcomplex<double> EigVec[1001]; // Eigenvectors stored sequentially.

  // Creating a complex matrix.

  nx = 10;
  n  =  nx*nx;
  CompMatrixA(nx, nnz, A, irow, pcol); 

  // Finding the four eigenvalues of A with largest magnitude
  // and the related eigenvectors.

  nconv = AREig(EigVal, EigVec, n, nnz, A, irow, pcol, 4);

  // Printing solution.

  Solution(nconv, n, nnz, A, irow, pcol, EigVal, EigVec);

} // main

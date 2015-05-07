/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE UCompReg.cc.
   Example program that illustrates how to solve a complex standard
   eigenvalue problem in regular mode using the ARluCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is obtained from the standard central difference
      discretization of the convection-diffusion operator
                    (Laplacian u) + rho*(du / dx)
      on the unit square [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      {nnz, irow, pcol, A}: matrix A data in CSC format.

   3) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      lcmatrxa.h       CompMatrixA, a function that generates matrix
                       A in CSC format.
      arunsmat.h       The ARumNonSymMatrix class definition.
      aruscomp.h       The ARluCompStdEig class definition.
      lcompsol.h       The Solution function.
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
#include "aruscomp.h"
#include "arunsmat.h"
#include "lcmatrxa.h"
#include "lcompsol.h"


int main()
{

  // Defining variables;

  int                nx;
  int                n;     // Dimension of the problem.
  int                nnz;   // Number of nonzero elements in A.
  int*               irow;  // pointer to an array that stores the row
                            // indices of the nonzeros in A.
  int*               pcol;  // pointer to an array of pointers to the
                            // beginning of each column of A in valA.
  arcomplex<double>* valA;  // pointer to an array that stores the
                            // nonzero elements of A.

  // Creating a complex matrix.

  nx = 10;
  n  = nx*nx;
  CompMatrixA(nx, nnz, valA, irow, pcol);
  ARumNonSymMatrix<arcomplex<double>, double> A(n, nnz, valA, irow, pcol);

  // Defining what we need: the four eigenvectors of A with largest magnitude.

  ARluCompStdEig<double> dprob(4L, A);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // main


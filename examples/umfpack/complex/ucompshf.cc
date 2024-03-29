/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE UCompShf.cc.
   Example program that illustrates how to solve a complex standard
   eigenvalue problem in shift and invert mode using the
   ARluCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and invert
      mode, where A is derived from the central difference discretization
      of the 1-dimensional convection-diffusion operator
                        (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      {nnz, irow, pcol, A}: matrix A data in CSC format.

   3) Library called by this example:

      The UMFPACK package is called by ARluCompStdEig to solve
      some linear systems involving (A-sigma*I). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      lcmatrxb.h       CompMatrixB, a function that generates matrix
                       A in CSC format.
      arunsmat.h       The ARumNonSymMatrix class definition.
      aruscomp.h       The ARluCompStdEig class definition.
      lcompsol.h       The Solution function.
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
#include "arunsmat.h"
#include "aruscomp.h"
#include "lcmatrxb.h"
#include "lcompsol.h"


int main()
{

  // Defining variables;

  int                n;     // Dimension of the problem.
  int                nnz;   // Number of nonzero elements in A.
  int*               irow;  // pointer to an array that stores the row
                            // indices of the nonzeros in A.
  int*               pcol;  // pointer to an array of pointers to the
                            // beginning of each column of A in valA.
  arcomplex<double>  rho;   // parameter used to generate A.
  arcomplex<double>* valA;  // pointer to an array that stores the
                            // nonzero elements of A.

  int nev = 4; // Number of requested eigenvalues.

  // Creating a complex matrix.

  n   = 100;
  rho = 10.0;
  CompMatrixB(n, rho, nnz, valA, irow, pcol);
  ARumNonSymMatrix<arcomplex<double>, double> A(n, nnz, valA, irow, pcol);

  // Defining what we need: the four eigenvectors of F nearest to 0.0.

  ARluCompStdEig<double> dprob(nev, A, arcomplex<double>(0.0, 0.0));

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

  int nconv = dprob.ConvergedEigenvalues();
  
  return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
} // main.


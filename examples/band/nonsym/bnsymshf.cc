/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE BNSymShf.cc.
   Example program that illustrates how to solve a real nonsymmetric
   standard eigenvalue problem in shift and invert mode using the
   ARluNonSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and invert
      mode, where A is derived from the standard central difference
      discretization of the 2-dimensional convection-diffusion operator
                       (Laplacian u) + rho*(du/dx)
      on a unit square with zero Dirichlet boundary conditions.
      The shift is a real number.

   2) Data structure used to represent matrix A:

      {ndiagL, ndiagU, A}: matrix A data in band format. The columns
      of A are stored sequentially in vector A. ndiagL and ndiagU
      supply the lower and upper bandwidth of A, respectively. 

   3) Library called by this example:

      The LAPACK package is called by ARluNonSymStdEig to solve
      some linear systems involving (A-sigma*I). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      --------------------------------------------
      bnmatrxa.h       BandMatrixA, a function that generates
                       matrix A in band format.
      arbnsmat.h       The ARbdNonSymMatrix class definition.
      arbsnsym.h       The ARluNonSymStdEig class definition.
      lnsymsol.h       The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "bnmatrxa.h"
#include "arbnsmat.h"
#include "arbsnsym.h"
#include "lnsymsol.h"


int main()
{

  // Defining variables;

  int     nx;
  int     n;       // Dimension of the problem.
  int     ndiagL;  // Lower bandwidth of A.
  int     ndiagU;  // Upper bandwidth of A.
  double  rho;     // Parameter used to define A.
  double* A;       // pointer to an array that stores the elements of A.

  int nev = 4; // Number of requested eigenvalues.

  // Creating a 100x100 matrix.

  nx  = 10;
  rho = 10.0;
  BandMatrixA(nx, rho, n, ndiagL, ndiagU, A);
  ARbdNonSymMatrix<double, double> matrix(n, ndiagL, ndiagU, A);

  // Defining what we need: the four eigenvectors of A nearest to 1000.0.

  ARluNonSymStdEig<double> dprob(nev, matrix, 1000.0, "LM", 20);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(matrix, dprob);

  int nconv = dprob.ConvergedEigenvalues();
  
  return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
} // main

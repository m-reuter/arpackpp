/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DSymGShf.cc.
   Example program that illustrates how to solve a real symmetric
   dense generalized eigenvalue problem in shift and invert mode 
   using the ARluSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A is the one dimensional discrete Laplacian 
      on the interval [0, 1], with zero Dirichlet boundary conditions, 
      and B is the mass matrix formed by using piecewise linear 
      elements on [0, 1].

   2) Data structure used to represent matrices A and B:

      Although A and B are very sparse in this example, they are 
      stored here as dense symmetric matrices. The lower triangular 
      part of A snd B is stored, by columns, in vectors A and B.

   3) Library called by this example:

      The LAPACK package is called by ARluSymGenEig to solve
      some linear systems involving (A-sigma*B).

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      dsmatrxb.h       DenseMatrixB, a function that generates
                       matrix A.
      dsmatrxc.h       DenseMatrixC, a function that generates 
                       matrix B.
      ardsmat.h        The ARdsSymMatrix class definition.
      ardgsym.h        The ARluSymGenEig class definition.
      lsymsol.h        The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "dsmatrxb.h"
#include "dsmatrxc.h"
#include "ardsmat.h"
#include "ardgsym.h"
#include "lsymsol.h"


int main()
{

  // Defining variables;

  int     n;       // Dimension of the problem.
  double* valA;    // pointer to an array that stores the elements of A.
  double* valB;    // pointer to an array that stores the elements of B.

  int nev = 4; // Number of requested eigenvalues.

  // Creating matrices A and B.

  n   = 100;
  DenseMatrixB(n, valA);
  ARdsSymMatrix<double> A(n, valA);

  DenseMatrixC(n, valB);
  ARdsSymMatrix<double> B(n, valB);

  // Defining what we need: the four eigenvectors nearest to 0.0.

  ARluSymGenEig<double> dprob('S', nev, A, B, 0.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

  int nconv = dprob.ConvergedEigenvalues();
  
  return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
} // main.


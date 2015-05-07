/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LSymShf.cc.
   Example program that illustrates how to solve a real symmetric
   standard eigenvalue problem in shift and invert mode using the
   ARluSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and  
      invert mode, where A is derived from the central difference 
      discretization of the one-dimensional Laplacian on [0, 1]
      with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      {nnz, irow, pcol, A}: lower triangular part of matrix A 
                            stored in CSC format.

   3) Library called by this example:

      The SuperLU package is called by ARluSymStdEig to solve
      some linear systems involving (A-sigma*I). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      --------------------------------------------
      lsmatrxa.h       SymmetricMatrixB, a function that generates
                       matrix A in CSC format.
      arlsmat.h        The ARluSymMatrix class definition.
      arlssym.h        The ARluSymStdEig class definition.
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

#include "lsmatrxb.h"
#include "arlsmat.h"
#include "arlssym.h"
#include "lsymsol.h"


int main()
{

  // Defining variables;

  int     n;          // Dimension of the problem.
  int     nnz;        // Number of nonzero elements in A.
  int*    irow;       // pointer to an array that stores the row
                      // indices of the nonzeros in A.
  int*    pcol;       // pointer to an array of pointers to the
                      // beginning of each column of A in vector A.
  double* A;          // pointer to an array that stores the
                      // nonzero elements of A.

  // Creating a 100x100 matrix.

  n = 100;
  SymmetricMatrixB(n, nnz, A, irow, pcol);
  ARluSymMatrix<double> matrix(n, nnz, A, irow, pcol);

  // Defining what we need: the four eigenvectors of A nearest to 1.0.

  ARluSymStdEig<double> dprob(4L, matrix, 1.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(matrix, dprob);

} // main

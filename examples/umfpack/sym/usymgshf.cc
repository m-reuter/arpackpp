/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE USymGShf.cc.
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in shift and invert mode using 
   the ARluSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are obtained from the finite element 
      discretization of the 1-dimensional discrete Laplacian
                                  d^2u / dx^2
      on the interval [0,1] with zero Dirichlet boundary conditions
      using piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: lower triangular part of matrix A 
                                  stored in CSC format.
      {nnzB, irowB, pcolB, valB}: lower triangular part of matrix B 
                                  stored in CSC format.

   3) Library called by this example:

      The UMFPACK package is called by ARluSymGenEig to solve
      some linear systems involving (A-sigma*B).

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lsmatrxc.h       SymmetricMatrixC, a function that generates
                       matrix A in CSC format.
      lsmatrxd.h       SymmetricMatrixD, a function that generates
                       matrix B in CSC format.
      arusmat.h        The ARumSymMatrix class definition.
      arugsym.h        The ARluSymGenEig class definition.
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

#include "lsmatrxc.h"
#include "lsmatrxd.h"
#include "arusmat.h"
#include "arugsym.h"
#include "lsymsol.h"


int main()
{

  int    n;              // Dimension of the problem.
  int    nnzA,   nnzB;   // Number of nonzero elements in A and B.
  int    *irowA, *irowB; // pointer to an array that stores the row
                         // indices of the nonzeros in A and B.
  int    *pcolA, *pcolB; // pointer to an array of pointers to the
                         // beginning of each column of A (B) in valA (valB).
  double *valA,  *valB;  // pointer to an array that stores the nonzero
                         // elements of A and B.

  // Creating matrices A and B.

  n = 100;
  SymmetricMatrixC(n, nnzA, valA, irowA, pcolA);
  ARumSymMatrix<double> A(n, nnzA, valA, irowA, pcolA);

  SymmetricMatrixD(n, nnzB, valB, irowB, pcolB);
  ARumSymMatrix<double> B(n, nnzB, valB, irowB, pcolB);

  // Defining what we need: the four eigenvectors nearest to 0.0.

  ARluSymGenEig<double> dprob('S', 4L, A, B, 0.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.


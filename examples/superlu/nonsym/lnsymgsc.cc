/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LNSymGSC.cc.
   Example program that illustrates how to solve a nonsymmetric
   generalized eigenvalue problem in complex shift and invert mode
   (taking the real part of OP*x) using the ARluNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in complex shift
      and inverse mode, where A is the tridiagonal matrix with 2 on the
      diagonal, -2 on the subdiagonal and 3 on the superdiagonal, and
      B is the tridiagonal matrix with 4 on the diagonal and 1 on the
      off-diagonals.
      The shift is a complex number.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.
      {nnzA, irowA, pcolA, valA}: matrix B data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by ARluNonSymGenEig to solve
      some complex linear systems involving (A-sigma*B). This is
      needed to implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      -----------------------------------------
      lnmatrxe.h       NonSymMatrixE, a function that generates
                       matrix A in CSC format.
      lnmatrxf.h       NonSymMatrixF, a function tha generates
                       matrix B in CSC format.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arlgnsym.h       The ARluNonSymGenEig class definition.
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

#include "lnmatrxe.h"
#include "lnmatrxf.h"
#include "arlnsmat.h"
#include "arlgnsym.h"
#include "lnsymsol.h"


int main()
{

  // Defining variables;

  int     n;               // Dimension of the problem.
  int     nnza,   nnzb;    // Number of nonzero elements in A and B.
  int     *irowa, *irowb;  // pointers to arrays that store the row
                           // indices of the nonzeros in A and B.
  int     *pcola, *pcolb;  // pointers to arrays of pointers to the
                           // beginning of each column of A and B in
                           // valA and valB.
  double  *valA,  *valB;   // pointers to arrays that store the
                           // nonzero elements of A and B.

  // Creating matrices A and B.

  n = 100;
  NonSymMatrixE(n, nnza, valA, irowa, pcola);
  ARluNonSymMatrix<double, double> A(n, nnza, valA, irowa, pcola);

  NonSymMatrixF(n, nnzb, valB, irowb, pcolb);
  ARluNonSymMatrix<double, double> B(n, nnzb, valB, irowb, pcolb);

  // Defining what we need: the four eigenvectors nearest to 0.4 + 0.6i.

  ARluNonSymGenEig<double> dprob(4L, A, B, 'R', 0.4, 0.6);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LNSymShf.cc.
   Example program that illustrates how to solve a real nonsymmetric
   standard eigenvalue problem in shift and invert mode using the
   ARluNonSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and invert
      mode, where A is derived from 2-D Brusselator Wave Model.
      The shift is a real number.

   2) Data structure used to represent matrix A:

      {nnz, irow, pcol, A}: matrix A data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by ARluNonSymStdEig to solve
      some linear systems involving (A-sigma*I). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      --------------------------------------------
      lnmatrxa.h       BrusselatorMatrix, a function that generates
                       matrix A in CSC format.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arlsnsym.h       The ARluNonSymStdEig class definition.
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

#include "lnmatrxa.h"
#include "arlnsmat.h"
#include "arlsnsym.h"
#include "lnsymsol.h"


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

  // Creating a 200x200 matrix.

  n = 200;
  BrusselatorMatrix(1.0, 0.004, 0.008, 2.0, 5.45, n, nnz, A, irow, pcol);
  ARluNonSymMatrix<double, double> BWM(n, nnz, A, irow, pcol);

  // Defining what we need: the four eigenvectors of BWM nearest to 0.0.

  ARluNonSymStdEig<double> dprob(4L, BWM, 0.0, "LM", 30L);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(BWM, dprob);

} // main.


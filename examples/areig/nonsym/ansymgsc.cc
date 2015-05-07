/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ANSymGSC.cc.
   Example program that illustrates how to solve a nonsymmetric
   generalized eigenvalue problem in complex shift and invert mode
   (taking the real part of OP*x) using the AREig function.

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

      The SuperLU package is called by AREig to solve some complex 
      linear systems involving (A-sigma*B). This is needed to 
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      -----------------------------------------
      lnmatrxe.h       NonSymMatrixE, a function that generates
                       matrix A in CSC format.
      lnmatrxf.h       NonSymMatrixF, a function tha generates
                       matrix B in CSC format.
      areig.h          The AREig function definition.
      ansymsol.h       The Solution function.

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
#include "areig.h"
#include "ansymsol.h"

int main()
{

  // Defining variables;

  int     n;               // Dimension of the problem.
  int     nconv;           // Number of "converged" eigenvalues.
  int     nnzA,   nnzB;    // Number of nonzero elements in A and B.
  int     *irowA, *irowB;  // pointers to arrays that store the row
                           // indices of the nonzeros in A and B.
  int     *pcolA, *pcolB;  // pointers to arrays of pointers to the
                           // beginning of each column of A and B in
                           // valA and ValB.
  float   *valA,  *valB;   // pointers to arrays that store the
                           // nonzero elements of A and B.
  float  EigValR[101];     // Real part of the eigenvalues.
  float  EigValI[101];     // Imaginary part of the eigenvalues.
  float  EigVec[1201];     // Eigenvectors stored sequentially.

  // Creating matrices A and B.

  n  = 100;  // Dimension of A and B.
  NonSymMatrixE(n, nnzA, valA, irowA, pcolA);
  NonSymMatrixF(n, nnzB, valB, irowB, pcolB);

  // Finding the four eigenvalues neares to 0.4 + 0.6I 
  // and the related eigenvectors.

  nconv = AREig(EigValR, EigValI, EigVec, n, nnzA, valA, 
                irowA, pcolA, nnzB, valB, irowB, pcolB, 
                'R', (float)0.4, (float)0.6, 4);

  // Printing solution.

  Solution(nconv, n, nnzA, valA, irowA, pcolA, nnzB,
           valB, irowB, pcolB, EigValR, EigValI, EigVec);

} // main.

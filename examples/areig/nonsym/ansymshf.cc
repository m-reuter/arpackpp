/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ANSymShf.cc
   Example program that illustrates how to solve a real nonsymmetric
   standard eigenvalue problem in shift and invert mode using the
   AREig function.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and invert
      mode, where A is derived from 2-D Brusselator Wave Model.
      The shift is a real number.

   2) Data structure used to represent matrix A:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by AREig to solve some linear 
      systems involving (A-sigma*I). This is needed to implement 
      the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      --------------------------------------------
      lnmatrxa.h       BrusselatorMatrix, a function that generates
                       matrix A in CSC format.
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

#include "lnmatrxa.h"
#include "areig.h"
#include "ansymsol.h"

int main()
{

  // Defining variables;

  int     n;            // Dimension of the problem.
  int     nconv;        // Number of "converged" eigenvalues.
  int     nnz;          // Number of nonzero elements in A.
  int*    irow;         // pointer to an array that stores the row
                        // indices of the nonzeros in A.
  int*    pcol;         // pointer to an array of pointers to the
                        // beginning of each column of A in vector A.
  double* A;            // pointer to an array that stores the
                        // nonzero elements of A.
  double  EigValR[101]; // Real part of the eigenvalues.
  double  EigValI[101]; // Imaginary part of the eigenvalues.
  double  EigVec[1201]; // Eigenvectors stored sequentially.

  // Creating a double precision matrix.

  n = 200;
  BrusselatorMatrix(1.0, 0.004, 0.008, 2.0, 5.45, n, nnz, A, irow, pcol);

  // Finding the four eigenvalues nearest to 0.0 and the
  // related eigenvectors.

  nconv = AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol, 0.0, 4);

  // Printing solution.

  Solution(nconv, n, nnz, A, irow, pcol, EigValR, EigValI, EigVec);

} // main.

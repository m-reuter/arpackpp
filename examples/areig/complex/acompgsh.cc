/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ACompGSh.cc.
   Example program that illustrates how to solve a complex
   generalized eigenvalue problem in shift and invert mode using
   the AREig function.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are derived from a finite element
      discretization of a 1-dimensional convection-diffusion operator
                         (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1], with zero boundary conditions, using
      piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.
      {nnzA, irowA, pcolA, valA}: matrix B data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by AREig to solve some linear 
      systems involving (A-sigma*B).

   4) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      lcmatrxe.h       CompMatrixE, a function that generates matrix
                       matrix A in CSC format.
      lcmatrxf.h       CompMatrixF, a function that generates matrix
                       B in CSC format.
      areig.h          The AREig function definition.
      acompsol.h       The Solution function.
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
#include "lcmatrxe.h"
#include "lcmatrxf.h"
#include "areig.h"
#include "acompsol.h"

int main()
{

  // Defining variables;

  int     n;                      // Dimension of the problem.
  int     nconv;                  // Number of "converged" eigenvalues.
  int     nnzA,   nnzB;           // Number of nonzero elements in A and B.
  int     *irowA, *irowB;         // pointers to arrays that store the row
                                  // indices of the nonzeros in A and B.
  int     *pcolA, *pcolB;         // pointers to arrays of pointers to the
                                  // beginning of each column of A and B in
                                  // valA and ValB.
  arcomplex<float> rho;           // parameter used in CompMatrixE.
  arcomplex<float> *valA, *valB;  // pointers to arrays that store the
                                  // nonzero elements of A and B.
  arcomplex<float> EigVal[101];   // Eigenvalues.
  arcomplex<float> EigVec[1001];  // Eigenvectors stored sequentially.

  // Creating matrices A and B.

  n   = 100;  // Dimension of A and B.
  rho = arcomplex<float>(10.0, 0.0);
  CompMatrixE(n, rho, nnzA, valA, irowA, pcolA);
  CompMatrixF(n, nnzB, valB, irowB, pcolB);

  // Finding the four eigenvalues nearest to 1.0 and the
  // related eigenvectors.

  nconv = AREig(EigVal, EigVec, n, nnzA, valA, irowA, pcolA, nnzB,
                valB, irowB, pcolB, arcomplex<float>(1.0, 0.0), 4);

  // Printing solution.

  Solution(nconv, n, nnzA, valA, irowA, pcolA, nnzB,
           valB, irowB, pcolB, EigVal, EigVec);

} // main.

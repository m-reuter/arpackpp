/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ANSymGRe.cc.
   Example program that illustrates how to solve a real 
   nonsymmetric generalized eigenvalue problem in regular 
   mode using the AREig function.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular
      mode, where A and B are derived from the finite element
      discretization of the 1-dimensional convection-diffusion operator
                        (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions
      using linear elements.

   2) Data structure used to represent matrices A and B:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.
      {nnzA, irowA, pcolA, valA}: matrix B data in CSC format.

   3) Library called by this example:

      The SuperLU package is called by AREig to solve some linear 
      systems involving B.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      lnmatrxc.h       StiffnessMatrix, a function that generates
                       matrix A in CSC format.
      lnmatrxd.h       MassMatrix, a function tha generates matrix
                       B in CSC format.
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

#include "lnmatrxc.h"
#include "lnmatrxd.h"
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
  double  rho;             // A parameter used in StiffnessMatrix.
  double  *valA,  *valB;   // pointers to arrays that store the
                           // nonzero elements of A and B.
  double EigValR[101];     // Real part of the eigenvalues.
  double EigValI[101];     // Imaginary part of the eigenvalues.
  double EigVec[1201];     // Eigenvectors stored sequentially.

  // Creating matrices A and B.

  n   = 100;  // Dimension of A and B.
  rho = 10.0;
  StiffnessMatrix(n, rho, nnzA, valA, irowA, pcolA);
  MassMatrix(n, nnzB, valB, irowB, pcolB);

  // Finding the four eigenvalues with largest magnitude 
  // and the related eigenvectors.

  nconv = AREig(EigValR, EigValI, EigVec, n, nnzA, valA, 
                irowA, pcolA, nnzB, valB, irowB, pcolB, 4); 

  // Printing solution.

  Solution(nconv, n, nnzA, valA, irowA, pcolA, nnzB,
           valB, irowB, pcolB, EigValR, EigValI, EigVec);

} // main


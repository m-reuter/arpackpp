/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE Simple.cc
   Simple example program that illustrates how to solve a real
   nonsymmetric standard eigenvalue problem in regular mode
   using the AREig function.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "lnmatrxc.h"
#include "areig.h"
#include <cmath>

int main()
{

  // Defining variables needed to store A in CSC format.

  int     n;     // Dimension of matrix.
  int     nnz;   // Number of nonzero elements in A.
  int*    irow;  // Row index of all nonzero elements of A.
  int*    pcol;  // Pointer to the beginning of each column (in irow and A).
  double* A;     // Nonzero elements of A.

  // Creating a double precision matrix.

  n = 100;
  StiffnessMatrix(n, 10.0, nnz, A, irow, pcol);

  // Defining AREig output variables.

  int     nconv;                       // Number of converged eigenvalues.
  double* EigValR = new double[201];   // Real part of the eigenvalues.
  double* EigValI = new double[201];   // Imaginary part of the eigenvalues.
  double* EigVec  = new double[1201];  // Eigenvectors.

  // Finding the five eigenvalues with largest magnitude
  // and the related eigenvectors.

  nconv = AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol, 5);

  // Printing eigenvalues.

  std::cout << "Eigenvalues:" << std::endl;
  for (int i=0; i<nconv; i++) {
    std::cout << "  lambda[" << (i+1) << "]: " << EigValR[i];
    if (EigValI[i]>=0.0) {
      std::cout << " + " << EigValI[i] << " I" << std::endl;
    }
    else {
      std::cout << " - " << fabs(EigValI[i]) << " I" << std::endl;
    }
  }
  std::cout << std::endl;

} // main

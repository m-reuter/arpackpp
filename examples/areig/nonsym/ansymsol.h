/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ANSymSol.cc
   Template functions that exemplify how to print information
   about nonsymmetric standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ANSYMSOL_H
#define ANSYMSOL_H

#include <cmath>
#include "blas1c.h"
#include "lapackc.h"
#include "arlnsmat.h"

template<class ARFLOAT, class ARINT>
void Solution(ARINT nconv, ARINT n, ARINT nnz, ARFLOAT A[], ARINT irow[], 
              ARINT pcol[], ARFLOAT EigValR[], ARFLOAT EigValI[], 
              ARFLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric eigen-problem
  on standard "std::cout" stream.
*/

{

  ARINT                              i;
  ARFLOAT*                           Ax;
  ARFLOAT*                           ResNorm;
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  std::cout << std::endl << std::endl << "Testing ARPACK++ function AREig" << std::endl;
  std::cout << "Real nonsymmetric eigenvalue problem: A*x - lambda*x \n \n";

  std::cout << "Dimension of the system            : " << n     << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv << std::endl << std::endl;

  // Printing eigenvalues.

  std::cout << "Eigenvalues:" << std::endl;

  for (i=0; i<nconv; i++) {
    std::cout << "  lambda[" << (i+1) << "]: " << EigValR[i];
    if (EigValI[i]>=0.0) {
      std::cout << " + " << EigValI[i] << " I" << std::endl;
    }
    else {
      std::cout << " - " << fabs(EigValI[i]) << " I" << std::endl;
    }
  }
  std::cout << std::endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Finding the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {

      if (EigValI[i]==0.0) { // Eigenvalue is real.

        matrix.MultMv(&EigVec[i*n], Ax);
        axpy(n, -EigValR[i], &EigVec[i*n], 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigValR[i]);

      }
      else {                 // Eigenvalue is complex.

        matrix.MultMv(&EigVec[i*n], Ax);
        axpy(n, -EigValR[i], &EigVec[i*n], 1, Ax, 1);
        axpy(n, EigValI[i], &EigVec[(i+1)*n], 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1);
        matrix.MultMv(&EigVec[(i+1)*n], Ax);
        axpy(n, -EigValI[i], &EigVec[i*n], 1, Ax, 1);
        axpy(n, -EigValR[i], &EigVec[(i+1)*n], 1, Ax, 1);
        ResNorm[i] = lapy2(ResNorm[i], nrm2(n, Ax, 1))/
                     lapy2(EigValR[i], EigValI[i]);
        ResNorm[i+1] = ResNorm[i];
        i++;

      }
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      std::cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << std::endl;
    }
    std::cout << std::endl;

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution.


template<class ARFLOAT, class ARINT>
void Solution(ARINT nconv, ARINT n, ARINT nnzA, ARFLOAT A[], ARINT irowA[],
              ARINT pcolA[], ARINT nnzB, ARFLOAT B[], ARINT irowB[], 
              ARINT pcolB[], ARFLOAT EigValR[], ARFLOAT EigValI[], 
              ARFLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric generalized
  eigen-problem on standard "std::cout" stream.
*/

{

  ARINT                              i;
  ARFLOAT                            *Ax;
  ARFLOAT                            *Bx, *Bx1;
  ARFLOAT                            *ResNorm;
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  std::cout << std::endl << std::endl << "Testing ARPACK++ function AREig" << std::endl;
  std::cout << "Real nonsymmetric generalized eigenvalue problem: A*x - lambda*B*x";
  std::cout << std::endl << std::endl;

  std::cout << "Dimension of the system            : " << n     << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv << std::endl << std::endl;

  // Printing eigenvalues.

  std::cout << "Eigenvalues:" << std::endl;

  for (i=0; i<nconv; i++) {
    std::cout << "  lambda[" << (i+1) << "]: " << EigValR[i];
    if (EigValI[i]>=0.0) {
      std::cout << " + " << EigValI[i] << " I" << std::endl;
    }
    else {
      std::cout << " - " << fabs(EigValI[i]) << " I" << std::endl;
    }
  }
  std::cout << std::endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    Bx      = new ARFLOAT[n];
    Bx1     = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {

      if (EigValI[i]==0.0) { // Eigenvalue is real.

        matrixA.MultMv(&EigVec[i*n], Ax);
        matrixB.MultMv(&EigVec[i*n], Bx);
        axpy(n, -EigValR[i], Bx, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigValR[i]);

      }
      else {                 // Eigenvalue is complex.

        matrixA.MultMv(&EigVec[i*n], Ax);
        matrixB.MultMv(&EigVec[i*n], Bx);
        matrixB.MultMv(&EigVec[(i+1)*n], Bx1);
        axpy(n, -EigValR[i], Bx, 1, Ax, 1);
        axpy(n, EigValI[i], Bx1, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1);
        matrixA.MultMv(&EigVec[(i+1)*n], Ax);
        axpy(n, -EigValI[i], Bx, 1, Ax, 1);
        axpy(n, -EigValR[i], Bx1, 1, Ax, 1);
        ResNorm[i] = lapy2(ResNorm[i], nrm2(n, Ax, 1))/
                     lapy2(EigValR[i], EigValI[i]);
        ResNorm[i+1] = ResNorm[i];
        i++;

      }
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << i << ") - lambda(" << i;
      std::cout << ")*B*x(" << i << ")||: " << ResNorm[i] << std::endl;
    }
    std::cout << std::endl;

    delete[] Ax;
    delete[] Bx;
    delete[] Bx1;
    delete[] ResNorm;

  }

} // Solution.


#endif // ANSYMSOL_H


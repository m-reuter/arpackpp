/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ASymSol.cc
   Template functions that exemplify how to print information
   about symmetric standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ASYMSOL_H
#define ASYMSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "arlsmat.h"

template<class ARFLOAT, class ARINT>
void Solution(ARINT nconv, ARINT n, ARINT nnz, ARFLOAT A[], ARINT irow[], 
              ARINT pcol[], char uplo, ARFLOAT EigVal[], ARFLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "cout" stream.
*/

{

  ARINT                  i;
  ARFLOAT*               Ax;
  ARFLOAT*               ResNorm;
  ARluSymMatrix<ARFLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric eigenvalue problem: A*x - lambda*x \n \n";

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Finding the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrix.MultMv(&EigVec[i*n], Ax);
      axpy(n, -EigVal[i], &EigVec[i*n], 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution.


template<class ARFLOAT, class ARINT>
void Solution(ARINT nconv, ARINT n, ARINT nnzA, ARFLOAT A[], ARINT irowA[],
              ARINT pcolA[], ARINT nnzB, ARFLOAT B[], ARINT irowB[], 
              ARINT pcolB[], char uplo, ARFLOAT EigVal[], ARFLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of symmetric generalized
  eigen-problem on standard "cout" stream.
*/

{

  ARINT                  i;
  ARFLOAT                *Ax, *Bx;
  ARFLOAT                *ResNorm;
  ARluSymMatrix<ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl << endl;

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    Bx      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrixA.MultMv(&EigVec[i*n], Ax);
      matrixB.MultMv(&EigVec[i*n], Bx);
      axpy(n, -EigVal[i], Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << i << ") - lambda(" << i;
      cout << ")*B*x(" << i << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution.


#endif // ASYMSOL_H


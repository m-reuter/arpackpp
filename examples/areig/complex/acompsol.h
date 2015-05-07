/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE ACompSol.cc
   Template functions that exemplify how to print information
   about complex standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ACOMPSOL_H
#define ACOMPSOL_H

#include <math.h>
#include "arcomp.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arlnsmat.h"

template<class ARFLOAT, class ARINT>
void Solution(ARINT nconv, ARINT n, ARINT nnz, arcomplex<ARFLOAT> A[], 
              ARINT irow[], ARINT pcol[], arcomplex<ARFLOAT> EigVal[], 
              arcomplex<ARFLOAT>* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric eigen-problems
  on standard "cout" stream.
*/

{

  ARINT                                         i;
  arcomplex<ARFLOAT>*                           Ax;
  ARFLOAT*                                      ResNorm;
  ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "complex standard eigenvalue problem: A*x - lambda*x \n \n";

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

    Ax      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrix.MultMv(&EigVec[i*n], Ax);
      axpy(n, -EigVal[i], &EigVec[i*n], 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/lapy2(real(EigVal[i]),imag(EigVal[i]));
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
void Solution(ARINT nconv, ARINT n, ARINT nnzA, arcomplex<ARFLOAT> A[], 
              ARINT irowA[], ARINT pcolA[], ARINT nnzB, 
              arcomplex<ARFLOAT> B[], ARINT irowB[], ARINT pcolB[],
              arcomplex<ARFLOAT> EigVal[], arcomplex<ARFLOAT>* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric generalized
  eigen-problem on standard "cout" stream.
*/

{

  ARINT                                        i;
  arcomplex<ARFLOAT>                           *Ax;
  arcomplex<ARFLOAT>                           *Bx;
  ARFLOAT                                      *ResNorm;
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Complex generalized eigenvalue problem: A*x - lambda*B*x \n \n";

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

    Ax      = new arcomplex<ARFLOAT>[n];
    Bx      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrixA.MultMv(&EigVec[i*n], Ax);
      matrixB.MultMv(&EigVec[i*n], Bx);
      axpy(n, -EigVal[i], Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/lapy2(real(EigVal[i]),imag(EigVal[i]));
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


#endif // ACOMPSOL_H


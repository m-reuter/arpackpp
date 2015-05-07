/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LSVDSol.h
   Template functions that exemplify how to print information
   about the singular value decomposition of a real nonsymmetric
   matrix.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LSVDSOL_H
#define LSVDSOL_H

#include <cmath>
#include "blas1c.h"
#include "lapackc.h"
#include "arssym.h"


template<class ARMATRIX, class ARFLOAT, class ARFOP>
void Solution(ARMATRIX &A, ARSymStdEig<ARFLOAT, ARFOP> &Prob)
/*
  Prints singular values and vectors on standard "std::cout" stream.
*/

{

  int     i, m, n, nAx, nconv;
  ARFLOAT *Ax;
  ARFLOAT *ResNorm;

  m     = A.nrows();
  n     = A.ncols();
  nAx   = (n>m)?n:m;
  nconv = Prob.ConvergedEigenvalues();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARSymStdEig \n";
  std::cout << "SVD problems: A = U*S*V'" << std::endl;

  std::cout << "Dimension of the system            : " << Prob.GetN()   << std::endl;
  std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev() << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv         << std::endl;
  std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv() << std::endl;
  std::cout << std::endl;

  if (Prob.EigenvaluesFound()) {

    // Printing singular values.

    std::cout << "Singular values:" << std::endl;
    for (i=0; i<nconv; i++) {
      std::cout << "  sigma[" << (i+1) << "]: " << Prob.Eigenvalue(i) << std::endl;
    }
    std::cout << std::endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*v - sigma*u ||
    // for the nconv accurately computed vectors u and v.

    Ax      = new ARFLOAT[nAx];
    ResNorm = new ARFLOAT[nconv+1];

    // Printing the residual norm || A*v - sigma*u ||
    // for the nconv accurately computed vectors u and v.

    for (i=0; i<nconv; i++) {

      A.MultMv(Prob.RawEigenvector(i)+m, Ax);
      axpy(m, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(m, Ax, 1)/fabs(Prob.Eigenvalue(i));

    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*v(" << (i+1) << ") - sigma(" << (i+1);
      std::cout << ")*u(" << (i+1) << ")||: " << ResNorm[i] << std::endl;
    }
    std::cout << std::endl;

    // Printing the residual norm || A'*u - sigma*v ||
    // for the nconv accurately computed vectors u and v.

    for (i=0; i<nconv; i++) {

      A.MultMtv(Prob.RawEigenvector(i), Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i)+m, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.Eigenvalue(i));

    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A'*u(" << (i+1) << ") - sigma(" << (i+1);
      std::cout << ")*v(" << (i+1) << ")||: " << ResNorm[i] << std::endl;
    }
    std::cout << std::endl;

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


#endif // LSVDSOL_H


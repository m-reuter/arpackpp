/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE CompGSol.h
   Template function that exemplify how to print information 
   about complex generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef COMPGSOL_H
#define COMPGSOL_H

#include "arcomp.h"
#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "argcomp.h"

template<class MATRA, class MATRB, class ARFOP, class ARFB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARCompGenEig<ARFLOAT, ARFOP, ARFB> &Prob)
/*
  Prints eigenvalues and eigenvectors of complex generalized
  eigen-problems on standard "cout" stream.
*/

{

  int                i, n, nconv, mode;
  ARFLOAT            *ResNorm;
  arcomplex<ARFLOAT> *Ax, *Bx;


  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARCompGenEig" << endl;
  cout << "Complex generalized eigenvalue problem: A*x - lambda*B*x" << endl;
  switch (mode) {
  case 2:
    cout << "Regular mode" << endl << endl;
    break;
  case 3:
    cout << "Shift and invert mode" << endl << endl;
  }

  cout << "Dimension of the system            : " << n              << endl;
  cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
  cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << endl;
  cout << "Number of iterations taken         : " << Prob.GetIter() << endl;
  cout << endl;

  if (Prob.EigenvaluesFound()) {

    // Printing eigenvalues.

    cout << "Eigenvalues:" << endl;
    for (i=0; i<nconv; i++) {
      cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << endl;
    }
    cout << endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new arcomplex<ARFLOAT>[n];
    Bx      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      B.MultMv(Prob.RawEigenvector(i),Bx);
      axpy(n, -Prob.Eigenvalue(i), Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/
                   lapy2(real(Prob.Eigenvalue(i)),imag(Prob.Eigenvalue(i)));
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*B*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    cout << endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution


#endif // COMPGSOL_H

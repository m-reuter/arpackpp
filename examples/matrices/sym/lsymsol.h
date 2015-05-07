/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LSymSol.h
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

#ifndef LSYMSOL_H
#define LSYMSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#ifdef ARLSMAT_H
#include "arlssym.h"
#include "arlgsym.h"
#elif defined ARUSMAT_H
#include "arussym.h"
#include "arugsym.h"
#elif defined ARDSMAT_H
#include "ardssym.h"
#include "ardgsym.h"
#else
#include "arbssym.h"
#include "arbgsym.h"
#endif


template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARluSymStdEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARluSymStdEig \n";
  cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << endl;
  switch (mode) {
  case 1:
    cout << "Regular mode" << endl;
    break;
  case 3:
    cout << "Shift and invert mode" << endl;
  }
  cout << endl;

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

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i), Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


template<class MATRA, class MATRB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARluSymGenEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric generalized
  eigen-problems on standard "cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax, *Bx, *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARluSymGenEig \n";
  cout << "Real symmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl;
  switch (mode) {
  case 2:
    cout << "Regular mode" << endl;
    break;
  case 3:
    cout << "Shift and invert mode" << endl;
    break;
  case 4:
    cout << "Buckling mode" << endl;
    break;
  case 5:
    cout << "Cayley mode" << endl;
  }
  cout << endl;

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

    Ax      = new ARFLOAT[n];
    Bx      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i), Ax);
      B.MultMv(Prob.RawEigenvector(i), Bx);
      axpy(n, -Prob.Eigenvalue(i), Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << i << ") - lambda(" << i;
      cout << ")*B*x(" << i << ")||: " << ResNorm[i] << "\n";
    }
    cout << "\n";

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution


#endif // LSYMSOL_H


/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LCompSol.h
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

#ifndef LCOMPSOL_H
#define LCOMPSOL_H

#include <math.h>
#include "arcomp.h"
#include "blas1c.h"
#include "lapackc.h"
#ifdef ARLNSMAT_H
#include "arlscomp.h"
#include "arlgcomp.h"
#elif defined ARUNSMAT_H
#include "aruscomp.h"
#include "arugcomp.h"
#elif defined ARDNSMAT_H
#include "ardscomp.h"
#include "ardgcomp.h"
#else
#include "arbscomp.h"
#include "arbgcomp.h"
#endif

template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARluCompStdEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues and eigenvectors of complex eigen-problems
  on standard "cout" stream.
*/

{

  int                i, n, nconv, mode;
  arcomplex<ARFLOAT> *Ax;
  ARFLOAT            *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARluCompStdEig \n";
  cout << "Complex eigenvalue problem: A*x - lambda*x" << endl;
  switch (mode) {
  case 1:
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

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/
                   lapy2(real(Prob.Eigenvalue(i)),imag(Prob.Eigenvalue(i)));
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << endl;
    }
    cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


template<class MATRA, class MATRB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARluCompGenEig<ARFLOAT> &Prob)
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

  cout << endl << endl; 
  cout << "Testing ARPACK++ class ARluCompGenEig \n" << endl;
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


#endif // LCOMPSOL_H


/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RSymGSol.h
   Printing eigenvalues of a symmetric generalized problem
   (ARrcSymGenEig version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef RSYMGSOL_H
#define RSYMGSOL_H

#include "arrgsym.h"

template<class ARFLOAT>
void Solution(ARrcSymGenEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues on standard "cout" stream.
*/

{

  int   i, n, nconv, mode;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARrcSymGenEig" << endl;
  cout << "Real symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
  switch (mode) {
  case 2:
    cout << "Regular mode" << endl << endl;
    break;
  case 3:
    cout << "Shift and invert mode" << endl << endl;
    break;
  case 4:
    cout << "Buckling mode" << endl << endl;
    break;
  case 5:
    cout << "Cayley mode" << endl << endl;
  }

  cout << "Dimension of the system            : " << n              << endl;
  cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
  cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << endl;
  cout << "Number of iterations taken         : " << Prob.GetIter() << endl;
  cout << endl;

  if (Prob.EigenvaluesFound()) {
    cout << "Eigenvalues:" << endl;
    for (i=0; i<nconv; i++) {
      cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << endl;
    }
    cout << endl;
  }

} // Solution


#endif // RSYMGSOL_H

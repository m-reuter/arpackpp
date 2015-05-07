/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RCompGSl.h
   Printing eigenvalues of a complex generalized problem
   (ARrcCompGenEig version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef RCOMPGSL_H
#define RCOMPGSL_H

#include "arrgcomp.h"

template<class ARFLOAT>
void Solution(ARrcCompGenEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues on standard "cout" stream.
*/

{

  int i, n, nconv, mode;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARrcCompGenEig" << endl;
  cout << "Complex eigenvalue problem: A*x - B*x*lambda" << endl;
  switch (mode) {
  case 2:
    cout << "Regular mode";
    break;
  case 3:
    cout << "Shift and invert mode";
  }
  cout << endl << endl;

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


#endif // RCOMPGSL_H

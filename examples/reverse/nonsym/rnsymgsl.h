/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RNSymGSl.h
   Printing eigenvalues of a nonsymmetric generalized problem
   (ARrcNonSymGenEig version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef RNSYMGSL_H
#define RNSYMGSL_H

#include "arrgnsym.h"

template<class ARFLOAT>
void Solution(ARrcNonSymGenEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues on standard "cout" stream.
*/

{

  int i, n, nconv, mode;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARrcNonSymGenEig" << endl;
  cout << "Real nonsymmetric generalized eigenvalue problem: A*x - B*x*lambda" << endl;
  switch (mode) {
  case 2:
    cout << "Regular mode";
    break;
  case 3:
    cout << "Shift and invert mode (using real part of OP)";
    break;
  case 4:
    cout << "Shift and invert mode (using imaginary part of OP)";
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


#endif // RNSYMGSL_H


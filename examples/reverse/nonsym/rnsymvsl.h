/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RNSymVSl.h
   Template functions that exemplify how to print information 
   about the singular value decomposition obtained using the
   ARrcNonSymStdEig function.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef RNSYMVSL_H
#define RNSYMVSL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "arrsnsym.h"

template<class ARFLOAT>
void Solution(ARrcNonSymStdEig<ARFLOAT> &Prob)
/*
  Prints singular values and singular vectors of nonsymmetric 
  real matrices on standard "cout" stream.
*/

{

  int  i, nconv;

  nconv = Prob.ConvergedEigenvalues();

  cout << endl << endl << "Testing ARPACK++ class ARrcNonSymStdEig \n";
  cout << "Singular value decomposition problem: (A'*A)*x - lambda*x" << endl;

  cout << "Dimension of the system              : " << Prob.GetN()   << endl;
  cout << "Number of 'requested' singular values: " << Prob.GetNev() << endl;
  cout << "Number of 'converged' singular values: " << nconv         << endl;
  cout << "Number of Arnoldi vectors generated  : " << Prob.GetNcv() << endl;
  cout << endl;

  if (Prob.EigenvaluesFound()) {

    // Printing singular values.

    cout << "Singular values:" << endl;
    for (i=0; i<nconv; i++) {
      cout << "  sigma[" << (i+1) << "]: ";
      cout << sqrt(Prob.EigenvalueReal(i)) << endl;
    }
    cout << endl;
  }

} // Solution


#endif // RNSYMVSL_H

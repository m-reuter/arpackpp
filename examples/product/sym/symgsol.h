/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE SymGSol.h
   Template functions that exemplify how to print information 
   about symmetric generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SYMGSOL_H
#define SYMGSOL_H

#include <math.h>
#include "blas1c.h"
#include "argsym.h"

template<class MATRA, class MATRB, class ARFOP, class ARFB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARSymGenEig<ARFLOAT, ARFOP, ARFB> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric generalized
  eigen-problems on standard "cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax, *Bx;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARSymGenEig \n";
  cout << "Real symmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl;
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
    ResNorm = new ARFLOAT[nconv];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      B.MultMv(Prob.RawEigenvector(i),Bx);
      axpy(n, -Prob.Eigenvalue(i), Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1) / fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*B*x(" << (i+1) << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution


#endif // SYMGSOL_H

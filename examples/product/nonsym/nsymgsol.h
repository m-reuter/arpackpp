/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymGSol.h
   Template functions that exemplify how to print information 
   about nonsymmetric generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NSYMGSOL_H
#define NSYMGSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "argnsym.h"

template<class MATRA, class MATRB, class ARFOP, class ARFB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARNonSymGenEig<ARFLOAT, ARFOP, ARFB> &Prob)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric generalized
  eigen-problems on standard "cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax;
  ARFLOAT *Bx, *Bx1;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARNonSymGenEig" << endl;
  cout << "Real nonsymmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl;
  switch (mode) {
  case 2:
    cout << "Regular mode" << endl;
    break;
  case 3:
    cout << "Shift and invert mode (using real part of OP)" << endl;
    break;
  case 4:
    cout << "Shift and invert mode (using imaginary part of OP)" << endl;
  }
  cout << endl;

  cout << "Dimension of the system            : " << n              << endl;
  cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
  cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << endl;
  cout << "Number of iterations taken         : " << Prob.GetIter() << endl;
  cout << endl;

  if (Prob.EigenvaluesFound()) {

    // Eigenvalues.

    cout << "Eigenvalues:" << endl;
    for (i=0; i<nconv; i++) {
      cout << "  lambda[" << (i+1) << "]: " << Prob.EigenvalueReal(i);
      if (Prob.EigenvalueImag(i)>=0.0) {
        cout << " + " << Prob.EigenvalueImag(i) << " I" << endl;
      }
      else {
        cout << " - " << fabs(Prob.EigenvalueImag(i)) << " I" << endl;
      }
    }
    cout << endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    Bx      = new ARFLOAT[n];
    Bx1     = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {

      if (Prob.EigenvalueImag(i)==0.0) { // Eigenvalue is real.

        A.MultMv(Prob.RawEigenvector(i), Ax);
        B.MultMv(Prob.RawEigenvector(i), Bx);
        axpy(n, -Prob.EigenvalueReal(i), Bx, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.EigenvalueReal(i));

      }
      else {                 // Eigenvalue is complex.

        A.MultMv(Prob.RawEigenvector(i), Ax);
        B.MultMv(Prob.RawEigenvector(i), Bx);
        B.MultMv(Prob.RawEigenvector(i+1), Bx1);
        axpy(n, -Prob.EigenvalueReal(i), Bx, 1, Ax, 1);
        axpy(n, Prob.EigenvalueImag(i), Bx1, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1);
        A.MultMv(Prob.RawEigenvector(i+1), Ax);
        axpy(n, -Prob.EigenvalueImag(i), Bx, 1, Ax, 1);
        axpy(n, -Prob.EigenvalueReal(i), Bx1, 1, Ax, 1);
        ResNorm[i] = lapy2(ResNorm[i], nrm2(n, Ax, 1))/
                     lapy2(Prob.EigenvalueReal(i), Prob.EigenvalueImag(i));
        ResNorm[i+1] = ResNorm[i];
        i++;

      }
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << i << ") - lambda(" << i;
      cout << ")*B*x(" << i << ")||: " << ResNorm[i] << "\n";
    }
    cout << "\n";

    delete[] Ax;
    delete[] Bx;
    delete[] Bx1;
    delete[] ResNorm;

  }

} // Solution


#endif // NSYMGSOL_H

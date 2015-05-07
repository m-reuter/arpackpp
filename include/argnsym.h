/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARGNSym.h.
   Arpack++ class ARNonSymGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARGNSYM_H
#define ARGNSYM_H

#include <stddef.h>
#include "arch.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arsnsym.h"
#include "argeig.h"
#include "arrgnsym.h"

template<class ARFLOAT, class ARFOP, class ARFB>
class ARNonSymGenEig:
  virtual public ARGenEig<ARFLOAT, ARFLOAT, ARFOP, ARFB>,
  virtual public ARNonSymStdEig<ARFLOAT, ARFOP>,
  virtual public ARrcNonSymGenEig<ARFLOAT>  {

 public:

 // a) Notation.

  typedef void (ARFB::* TypeBx)(ARFLOAT[], ARFLOAT[]);


 protected:

 // b) Protected variables:

  ARFB    *objA;      // Object that has MultAx as a member function.
  TypeBx  MultAx;     // Function that evaluates the product A*x.


 // c) Protected functions:

  void RecoverEigenvalues();
  // Uses Rayleigh quotient to recover eigenvalues of the original
  // problem when shift is complex.

  virtual void Copy(const ARNonSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // d) Public functions:

 // d.1) Functions that allow changes in problem parameters.

  virtual void SetShiftInvertMode(ARFLOAT sigmaRp, ARFOP* objOPp,
                                  void (ARFOP::* MultOPxp)(ARFLOAT[],ARFLOAT[]));
  // Turns the problem to real shift-and-invert mode with sigmaRp as shift.

  virtual void SetComplexShiftMode(char partp, ARFLOAT sigmaRp, 
                                   ARFLOAT sigmaIp, ARFOP* objOPp, 
                                   void (ARFOP::* MultOPxp)(ARFLOAT[],ARFLOAT[]), 
                                   ARFB* objAp,
                                   void (ARFB::* MultAxp)(ARFLOAT[],ARFLOAT[]));
  // Turns the problem to complex shift-and-invert mode with shift
  // defined by sigmaRp and sigmaIp. MultAx is used to obtain eigenvalues.


 // d.2) Functions that perform all calculations in one step.

  virtual int FindEigenvalues();
  // Determines nev approximated eigenvalues of the given eigen-problem.

  virtual int FindEigenvectors(bool schurp = false);
  // Determines nev approximated eigenvectors of the given eigen-problem
  // Optionally also determines nev Schur vectors that span the desired
  // invariant subspace.

  virtual int FindSchurVectors();
  // Determines nev Schur vectors that span the desired invariant subspace.
  // Redefined in ARSymEig.


 // d.3) Constructors and destructor.

  ARNonSymGenEig() { part = 'R'; }
  // Short constructor (Does nothing but calling base classes constructors).

  ARNonSymGenEig(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
                 ARFB* objBp, void (ARFB::* MultBxp)(ARFLOAT[], ARFLOAT[]),
                 char* whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
                 int maxitp = 0, ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARNonSymGenEig(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
                 ARFB* objBp, void (ARFB::* MultBxp)(ARFLOAT[], ARFLOAT[]),
                 ARFLOAT sigmap, char* whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0, ARFLOAT* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (real shift and invert mode).

  ARNonSymGenEig(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]), ARFB* objAp,
                 void (ARFB::* MultAxp)(ARFLOAT[], ARFLOAT[]), ARFB* objBp,
                 void (ARFB::* MultBxp)(ARFLOAT[], ARFLOAT[]), char partp,
                 ARFLOAT sigmaRp, ARFLOAT sigmaIp, char* whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (complex shift and invert mode).

  ARNonSymGenEig(const ARNonSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARNonSymGenEig() { }
  // Destructor.

 // e) Operators.

  ARNonSymGenEig& operator=(const ARNonSymGenEig& other);
  // Assignment operator.

}; // class ARNonSymGenEig.


// ------------------------------------------------------------------------ //
// ARNonSymGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARFOP, class ARFB>
inline void ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::
Copy(const ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>& other)
{

  ARGenEig<ARFLOAT, ARFLOAT, ARFOP, ARFB>::Copy(other);
  objA   = other.objA;
  MultAx = other.MultAx;
  part   = other.part;

} // Copy.


template<class ARFLOAT, class ARFOP, class ARFB>
void ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::RecoverEigenvalues()
{

  int    j, ColJ, ColJp1;
  ARFLOAT  numr, numi, denr, deni;
  ARFLOAT* Ax;

  Ax = new ARFLOAT[n];

  for (j=0; j<nconv; j++) {

    ColJ   = j*n;
    ColJp1 = ColJ+n;

    if (EigValI[j] == (ARFLOAT)0.0) {

      // Eigenvalue is real. Computing EigVal = x'(Ax)/x'(Mx).

      (objB->*MultAx)(&EigVec[ColJ], Ax);
      numr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      (objB->*MultBx)(&EigVec[ColJ], Ax);
      denr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      EigValR[j] =  numr / denr;

    }
    else {

      // Eigenvalue is complex.

      // Computing x'(Ax).

      (objB->*MultAx)(&EigVec[ColJ], Ax);
      numr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      numi = dot(n, &EigVec[ColJp1], 1, Ax, 1);
      (objB->*MultAx)(&EigVec[ColJp1], Ax);
      numr = numr + dot(n, &EigVec[ColJp1], 1, Ax, 1);
      numi = -numi + dot(n, &EigVec[ColJ], 1, Ax, 1);

      // Computing x'(Mx).

      (objB->*MultBx)(&EigVec[ColJ], Ax);
      denr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      deni = dot(n, &EigVec[ColJp1], 1, Ax, 1);
      (objB->*MultBx)(&EigVec[ColJp1], Ax);
      denr = denr + dot(n, &EigVec[ColJp1], 1, Ax, 1);
      deni = -deni + dot(n, &EigVec[ColJ], 1, Ax, 1);

      // Computing the first eigenvalue of the conjugate pair.

      EigValR[j] = (numr*denr+numi*deni) / lapy2(denr, deni);
      EigValI[j] = (numi*denr-numr*deni) / lapy2(denr, deni);

      // Getting the second eigenvalue of the conjugate pair by taking
      // the conjugate of the first.

      EigValR[j+1] = EigValR[j];
      EigValI[j+1] = -EigValI[j];
      j++;

    }

  }

  delete[] Ax;

} // RecoverEigenvalues.


template<class ARFLOAT, class ARFOP, class ARFB>
inline void ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::
SetShiftInvertMode(ARFLOAT sigmaRp, ARFOP* objOPp,
                   void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]))
{

  part    = 'R';
  objOP   = objOPp;
  MultOPx = MultOPxp;
  ChangeShift(sigmaRp);

} // SetShiftInvertMode.


template<class ARFLOAT, class ARFOP, class ARFB>
inline void ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::
SetComplexShiftMode(char partp, ARFLOAT sigmaRp, ARFLOAT sigmaIp, 
                    ARFOP* objOPp, 
                    void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
                    ARFB* objAp, void (ARFB::* MultAxp)(ARFLOAT[], ARFLOAT[]))
{

  objOP   = objOPp;
  MultOPx = MultOPxp;
  objA    = objAp;
  MultAx  = MultAxp;
  part    = CheckPart(partp);
  ChangeShift(sigmaRp, sigmaIp);

} // SetComplexShiftMode.


template<class ARFLOAT, class ARFOP, class ARFB>
inline int ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::FindEigenvalues()
{

  nconv = ARStdEig<ARFLOAT, ARFLOAT, ARFOP>::FindEigenvalues();
  if (sigmaI != 0.0) RecoverEigenvalues();
  return nconv;

} // FindEigenvalues.


template<class ARFLOAT, class ARFOP, class ARFB>
inline int ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::FindEigenvectors(bool schurp)
{

  nconv = ARStdEig<ARFLOAT, ARFLOAT, ARFOP>::FindEigenvectors(schurp);
  if (sigmaI != 0.0) RecoverEigenvalues();
  return nconv;

} // FindEigenvectors.


template<class ARFLOAT, class ARFOP, class ARFB>
int ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::FindSchurVectors()
{

  nconv = ARStdEig<ARFLOAT, ARFLOAT, ARFOP>::FindSchurVectors();
  if (sigmaI != 0.0) RecoverEigenvalues();
  return nconv;

} // FindSchurVectors.


template<class ARFLOAT, class ARFOP, class ARFB>
inline ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::
ARNonSymGenEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
               ARFB* objBp, void (ARFB::* MultBxp)(ARFLOAT[], ARFLOAT[]),
               char* whichp, int ncvp, ARFLOAT tolp, int maxitp,
               ARFLOAT* residp, bool ishiftp)

{

  part = 'R';                // Considering mode = 3 in ChangeShift.
  NoShift();
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT, class ARFOP, class ARFB>
inline ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::
ARNonSymGenEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
               ARFB* objBp, void (ARFB::* MultBxp)(ARFLOAT[], ARFLOAT[]),
               ARFLOAT sigmap, char* whichp, int ncvp,
               ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  SetShiftInvertMode(sigmap, objOPp, MultOPxp);
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);


} // Long constructor (real shift and invert mode).


template<class ARFLOAT, class ARFOP, class ARFB>
inline ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::
ARNonSymGenEig(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
               ARFB* objAp, void (ARFB::* MultAxp)(ARFLOAT[], ARFLOAT[]),
               ARFB* objBp, void (ARFB::* MultBxp)(ARFLOAT[], ARFLOAT[]),
               char partp, ARFLOAT sigmaRp, ARFLOAT sigmaIp,
               char* whichp, int ncvp, ARFLOAT tolp, int maxitp,
               ARFLOAT* residp, bool ishiftp)

{

  SetComplexShiftMode(partp, sigmaRp, sigmaIp, objOPp,
                      MultOPxp, objAp, MultAxp);
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT, class ARFOP, class ARFB>
ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>& ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>::
operator=(const ARNonSymGenEig<ARFLOAT, ARFOP, ARFB>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARGNSYM_H


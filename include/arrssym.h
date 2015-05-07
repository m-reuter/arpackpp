/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARRSSym.h.
   Arpack++ class ARrcSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRSSYM_H
#define ARRSSYM_H

#include <stddef.h>
#include "arch.h"
#include "arerror.h"
#include "debug.h"
#include "arrseig.h"
#include "saupp.h"
#include "seupp.h"


template<class ARFLOAT>
class ARrcSymStdEig: public virtual ARrcStdEig<ARFLOAT, ARFLOAT> {

 protected:

 // a) Protected functions:

 // a.1) Memory control functions.

  void WorkspaceAllocate();
  // Allocates workspace for symmetric problems.


 // a.2) Functions that handle original FORTRAN ARPACK code.

  void Aupp();
  // Interface to FORTRAN subroutines SSAUPD and DSAUPD.

  void Eupp();
  // Interface to FORTRAN subroutines SSEUPD and DSEUPD.


 // a.3) Functions that check user defined parameters.

  char* CheckWhich(char* whichp);
  // Determines if the value of variable "which" is valid.


 public:

 // b) Public functions:

 // b.1) Trace functions.

  void Trace(const int digit = -5, const int getv0 = 0, const int aupd = 1,
             const int aup2 = 0,  const int aitr = 0,  const int eigt = 0,
             const int apps = 0,  const int gets = 0,  const int eupd = 0) 
  { 
    sTraceOn(digit, getv0, aupd, aup2, aitr, eigt, apps, gets, eupd);
  }
  // Turns on trace mode. 


 // b.2) Functions that permit step by step execution of ARPACK.

  ARFLOAT* PutVector();
  // When ido = -1, 1 or 2 and the user must perform a product in the form
  // y <- M*x, this function indicates where to store y. When ido = 3, this
  // function indicates where to store the shifts.


 // b.3) Functions that perform all calculations in one step.

  int FindSchurVectors() {
    throw ArpackError(ArpackError::SCHUR_UNDEFINED, "FindSchurVectors");
    return 0;  // Only to avoid warning messages emitted by some compilers.
  }
  // For symmetric problems, Schur vectors are eigenvectors.

  int Eigenvalues(ARFLOAT* &EigValp, bool ivec = false, bool ischur = false);
  // Overrides array EigValp with the eigenvalues of the problem.
  // Also calculates eigenvectors and Schur vectors if requested.

  int EigenValVectors(ARFLOAT* &EigVecp, ARFLOAT* &EigValp, 
                      bool ischur = false);
  // Overrides array EigVecp sequentially with the eigenvectors of the
  // given eigen-problem. Also stores the eigenvalues in EigValp.
  // Calculates Schur vectors if requested.


 // b.4) Functions that return elements of vectors and matrices.

  ARFLOAT Eigenvalue(int i);
  // Provides i-eth eigenvalue.

  ARFLOAT Eigenvector(int i, int j);
  // Provides element j of the i-eth eigenvector.


 // b.5) Functions that use STL vector class.

#ifdef STL_VECTOR_H

  vector<ARFLOAT>* StlEigenvalues(bool ivec = false, bool ischur = false);
  // Calculates the eigenvalues and stores them in a single STL vector.
  // Also calculates eigenvectors and Schur vectors if requested.

  vector<ARFLOAT>* StlEigenvector(int i);
  // Returns the i-th eigenvector in a STL vector.

#endif // #ifdef STL_VECTOR_H.


 // b.6) Constructors and destructor.

  ARrcSymStdEig() { }
  // Short constructor.

  ARrcSymStdEig(int np, int nevp, char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0, ARFLOAT* residp = NULL,
                bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcSymStdEig(int np, int nevp, ARFLOAT sigmap, char* whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcSymStdEig(const ARrcSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcSymStdEig() { }
  // Destructor.

 // c) Operators.

  ARrcSymStdEig& operator=(const ARrcSymStdEig& other);
  // Assignment operator.

}; // class ARrcSymStdEig.


// ------------------------------------------------------------------------ //
// ARrcSymStdEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARrcSymStdEig<ARFLOAT>::WorkspaceAllocate()
{

  lworkl  = ncv*(ncv+9);
  lworkv  = 0;
  lrwork  = 0;
  workl   = new ARFLOAT[lworkl+1];

} // WorkspaceAllocate.


template<class ARFLOAT>
inline void ARrcSymStdEig<ARFLOAT>::Aupp()
{

  saupp(ido, bmat, n, which, nev, tol, resid, ncv, V, n,
        iparam, ipntr, workd, workl, lworkl, info);

} // Aupp.


template<class ARFLOAT>
inline void ARrcSymStdEig<ARFLOAT>::Eupp()
{

  seupp(rvec, HowMny, EigValR, EigVec, n, sigmaR, bmat,
        n, which, nev, tol, resid, ncv, V, n, iparam,
        ipntr, workd, workl, lworkl, info);

} // Eupp.


template<class ARFLOAT>
char* ARrcSymStdEig<ARFLOAT>::CheckWhich(char* whichp)
{

  switch (whichp[0]) {
  case 'B':                       // The options are: BE, ...
    return "BE";
  case 'L':                       // LA, LM, ...
  case 'S':                       // SA, SM.
    switch (whichp[1]){
    case 'A':
    case 'M':
      return whichp;
    }
  default:
    throw ArpackError(ArpackError::WHICH_UNDEFINED);
  }

} // CheckWhich.


template<class ARFLOAT>
ARFLOAT* ARrcSymStdEig<ARFLOAT>::PutVector()
{

  switch (ido) {
  case -1:
  case  1:                    // Returning OP*x.
  case  2:
    return &workd[ipntr[2]];  // Returning B*x.
  case  3:
    return &workl[ipntr[11]]; // Returning shifts.
  default:
    throw ArpackError(ArpackError::CANNOT_PUT_VECTOR, "PutVector");
  }

} // PutVector.


template<class ARFLOAT>
int ARrcSymStdEig<ARFLOAT>::
Eigenvalues(ARFLOAT* &EigValp, bool ivec, bool ischur)
{

  if (ValuesOK) {                      // Eigenvalues are available.
    if (EigValp == NULL) {             // Moving eigenvalues.
      EigValp  = EigValR;
      EigValR  = NULL;
      newVal   = false;
      ValuesOK = false;
    }
    else {                             // Copying eigenvalues.
      copy(nconv,EigValR,1,EigValp,1);
    }
  }
  else {                               // Eigenvalues are not available.
    if (newVal) {
      delete[] EigValR;
      newVal = false;
    }
    if (EigValp == NULL) {
      try { EigValp = new ARFLOAT[ValSize()]; }
      catch (ArpackError) { return 0; }
    }
    EigValR = EigValp;
    if (ivec) {                        // Finding eigenvalues and eigenvectors.
      nconv = FindEigenvectors(ischur);
    }
    else {                             // Finding eigenvalues only.
      nconv = FindEigenvalues();
    }
    EigValR = NULL;
  }
  return nconv;

} // Eigenvalues(EigValp, ivec, ischur).


template<class ARFLOAT>
int ARrcSymStdEig<ARFLOAT>::
EigenValVectors(ARFLOAT* &EigVecp, ARFLOAT* &EigValp, bool ischur)
{

  if (ValuesOK) {                  // Eigenvalues are already available .
    nconv = Eigenvalues(EigValp, false);
    nconv = Eigenvectors(EigVecp, ischur);
  }
  else {                           // Eigenvalues and vectors are not available.
    try {
      if (EigVecp == NULL) EigVecp = new ARFLOAT[ValSize()*n];
      if (EigValp == NULL) EigValp = new ARFLOAT[ValSize()];
    }
    catch (ArpackError) { return 0; }
    if (newVec) {
      delete[] EigVec;
      newVec = false;
    }
    if (newVal) {
      delete[] EigValR;
      newVal = false;
    }
    EigVec  = EigVecp;
    EigValR = EigValp;
    nconv   = FindEigenvectors(ischur);
    EigVec  = NULL;
    EigValR = NULL;
  }
  return nconv;

} // EigenValVectors(EigVecp, EigValp, ischur).


template<class ARFLOAT>
inline ARFLOAT ARrcSymStdEig<ARFLOAT>::Eigenvalue(int i)
{

  // Returning i-eth eigenvalue.

  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "Eigenvalue(i)");
  }
  else if ((i>=nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvalue(i)");
  }
  return EigValR[i];

} // Eigenvalue(i).


template<class ARFLOAT>
inline ARFLOAT ARrcSymStdEig<ARFLOAT>::Eigenvector(int i, int j)
{

  // Returning element j of i-eth eigenvector.

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "Eigenvector(i,j)");
  }
  else if ((i>=nconv)||(i<0)||(j>=n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvector(i,j)");
  }
  return EigVec[i*n+j];

} // Eigenvector(i,j).


#ifdef STL_VECTOR_H

template<class ARFLOAT>
inline vector<ARFLOAT>* ARrcSymStdEig<ARFLOAT>::
StlEigenvalues(bool ivec, bool ischur)
{

  // Returning the eigenvalues in a STL vector.

  vector<ARFLOAT>* StlEigValR;
  ARFLOAT*         ValPtr;

  try { StlEigValR = new vector<ARFLOAT>(ValSize()); }
  catch (ArpackError) { return NULL; }
  ValPtr = StlEigValR->begin();
  nconv = Eigenvalues(ValPtr, ivec, ischur);
  return StlEigValR;

} // StlEigenvalues.


template<class ARFLOAT>
inline vector<ARFLOAT>* ARrcSymStdEig<ARFLOAT>::StlEigenvector(int i)
{

  // Returning the i-th eigenvector in a STL vector.

  vector<ARFLOAT>* Vec;

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvector(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvector(i)");
  }
  try {
    Vec = new vector<ARFLOAT>(&EigVec[i*n], &EigVec[(i+1)*n]);
  }
  catch (ArpackError) { return NULL; }
  return Vec;

} // StlEigenvector(i).

#endif // #ifdef STL_VECTOR_H.


template<class ARFLOAT>
inline ARrcSymStdEig<ARFLOAT>::
ARrcSymStdEig(int np, int nevp, char* whichp, int ncvp,
              ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARrcSymStdEig<ARFLOAT>::
ARrcSymStdEig(int np, int nevp, ARFLOAT sigmap, char* whichp,
              int ncvp, ARFLOAT tolp, int maxitp, ARFLOAT* residp,
              bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARrcSymStdEig<ARFLOAT>& ARrcSymStdEig<ARFLOAT>::
operator=(const ARrcSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRSSYM_H


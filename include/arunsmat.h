/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUNSMat.h.
   Arpack++ class ARumNonSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arunspen.h"

#ifndef ARUNSMAT_H
#define ARUNSMAT_H

#include <algorithm>
#include <cstddef>
#include <string>
#include "arch.h"
#include "armat.h"
#include "arspmat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "umfpackc.h"

template<class ARTYPE, class ARFLOAT> class ARumNonSymPencil;

template<class ARTYPE, class ARFLOAT>
class ARumNonSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARumNonSymPencil<ARTYPE, ARFLOAT>;
  friend class ARumNonSymPencil<ARFLOAT, ARFLOAT>;

 protected:

  double  control[UMFPACK_CONTROL];
  double  info[UMFPACK_INFO];
  void*   Numeric;
  bool    factored;

  // The input matrix
  ARSparseMatrix<ARTYPE>* mat;

  virtual void Copy(const ARumNonSymMatrix& other);

  void ClearMem();

  void SubtractAsI(ARTYPE sigma);

  void Check(int status);

 private:

  // Internal reference to current matrix (either input matrix or AsI)
  ARSparseMatrix<ARTYPE>* pA;

  // Internal matrix storing A - s I
  ARSparseMatrix<ARTYPE>* AsI;

 public:

  int nzeros() { return mat->nzeros(); }

  bool IsSymmetric() { return false /*bool(icntl[5])*/; }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultMtv(ARTYPE* v, ARTYPE* w);

  void MultMtMv(ARTYPE* v, ARTYPE* w);

  void MultMMtv(ARTYPE* v, ARTYPE* w);

  void Mult0MMt0v(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
                    double thresholdp = 0.1, bool check = true,
                    bool owner = false); // Square.

  void DefineMatrix(int mp, int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
                    bool check = true, bool owner = false); // Rectangular.

  ARumNonSymMatrix(): ARMatrix<ARTYPE>(), factored(false), Numeric(nullptr), mat(nullptr), AsI(nullptr)
  {
  }
  // Short constructor that does nothing.

  ARumNonSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
                   double thresholdp = 0.1, bool check = true);
  // Long constructor (square matrix).

  ARumNonSymMatrix(int mp, int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp);
  // Long constructor (rectangular matrix).

  ARumNonSymMatrix(const std::string& name, double thresholdp = 0.1, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARumNonSymMatrix(const ARumNonSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumNonSymMatrix() { ClearMem(); }
  // Destructor.

  ARumNonSymMatrix& operator=(const ARumNonSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumNonSymMatrix member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::ClearMem()
{

  if (factored && Numeric)
  {
    umfpack_free_numeric<ARTYPE>(&Numeric);
    Numeric = nullptr;
  }

  if (mat) { delete mat; mat = nullptr; }
  if (AsI) { delete AsI; AsI = nullptr; }

  pA = nullptr;

} // ClearMem.


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::
Copy(const ARumNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  this->m         = other.m;
  this->n         = other.n;
  this->defined   = other.defined;
  factored  = other.factored;

  mat->Copy(*other.mat);

  // Returning from here if "other" was not initialized.

  if (!this->defined) return;

  // Copying arrays with static dimension.

  for (int i = 0; i < UMFPACK_CONTROL; i++) control[i] = other.control[i];
  for (int i = 0; i < UMFPACK_INFO; i++) info[i] = other.info[i];

  // Returning from here if "other" was not factored.

  if (!factored) return;

  factored = false;
  
} // Copy.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::SubtractAsI(ARTYPE sigma)
{
  if (!AsI)
  {
      int ndiag = mat->DiagIndices();
      int nz = mat->nzeros();

      AsI = new ARSparseMatrix<ARTYPE>(this->m, this->n, nz + this->n - ndiag);
  }

  AsI->Copy(*mat);

  if (sigma != (ARTYPE)0)
  {
    AsI->AddDiag(-sigma);
  }

  pA = AsI;

}


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::Check(int status)
{

  // status = info[0]
  if (status == UMFPACK_ERROR_out_of_memory)  {
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARumNonSymMatrix::FactorA");
  }
  else if (status == UMFPACK_ERROR_invalid_matrix)  {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumNonSymMatrix::FactorA");
  }
  else if (status == UMFPACK_WARNING_singular_matrix) {
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARumNonSymMatrix::FactorA");
  }
  else if (status != UMFPACK_OK) {
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARumNonSymMatrix::FactorA");
  }

} // Check.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,"ARumNonSymMatrix::FactorA");
  }

  // Quitting the function if A is not square.

  if (this->m != this->n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymMatrix::FactorA");
  }

  // Decomposing A.

  void *Symbolic;

  auto ap = mat->pcol();
  auto ai = mat->irow();
  auto ax = mat->values();

  Check(umfpack_symbolic(this->m, this->n, ap, ai, ax, &Symbolic, control, info));
  Check(umfpack_numeric(ap, ai, ax, Symbolic, &Numeric, control, info));

  umfpack_free_symbolic<ARTYPE>(&Symbolic);

  factored = true;

  pA = mat;

} // FactorA.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymMatrix::FactorAsI");
  }

  // Quitting the function if A is not square.

  if (this->m != this->n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymMatrix::FactorAsI");
  }

  // Subtracting sigma*I from A (this will allocate AsI).

  SubtractAsI(sigma);

  // Decomposing AsI.

  void *Symbolic;

  auto ap = AsI->pcol();
  auto ai = AsI->irow();
  auto ax = AsI->values();

  Check(umfpack_symbolic(this->m, this->n, ap, ai, ax, &Symbolic, control, info));
  Check(umfpack_numeric(ap, ai, ax, Symbolic, &Numeric, control, info));

  umfpack_free_symbolic<ARTYPE>(&Symbolic);

  factored = true;

  pA = AsI;

} // FactorAsI.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::MultMv(ARTYPE* v, ARTYPE* w)
{

  int    i,j;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumNonSymMatrix::MultMv");
  }

  auto ax = pA->values();
  auto ap = pA->pcol();
  auto ai = pA->irow();

  // Determining w = M.v.

  for (i = 0; i != this->m; i++) w[i]=(ARTYPE)0;

  for (i = 0; i != this->n; i++) {
    t = v[i];
    for (j=ap[i]; j!=ap[i+1]; j++) {
      w[ai[j]] += t*ax[j];
    }
  }

} // MultMv.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::MultMtv(ARTYPE* v, ARTYPE* w)
{

  int    i,j;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,"ARumNonSymMatrix::MultMtv");
  }

  // Determining w = M'.v.

  auto ax = pA->values();
  auto ap = pA->pcol();
  auto ai = pA->irow();

  for (i = 0; i != this->n; i++) {
    t = (ARTYPE)0;
    for (j = ap[i]; j != ap[i+1]; j++) {
      t += v[ai[j]]*ax[j];
    }
    w[i] = t;
  }

} // MultMtv.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::MultMtMv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE* t = new ARTYPE[this->m];

  MultMv(v,t);
  MultMtv(t,w);

  delete[] t;

} // MultMtMv.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::MultMMtv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE* t = new ARTYPE[this->n];

  MultMtv(v,t);
  MultMv(t,w);

  delete[] t;

} // MultMMtv.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::Mult0MMt0v(ARTYPE* v, ARTYPE* w)
{

  MultMv(&v[this->m],w);
  MultMtv(v,&w[this->m]);

} // Mult0MMt0v.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARumNonSymMatrix::MultInvv");
  }

  auto ap = pA->pcol();
  auto ai = pA->irow();
  auto ax = pA->values();

  // Solving A.w = v (or AsI.w = v).

  int status = umfpack_solve(UMFPACK_A, ap, ai, ax, w, v, Numeric, control, info);

  if (status != UMFPACK_OK)
      throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumNonSymMatrix::MultInvv");

} // MultInvv.


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
             double thresholdp, bool check, bool owner)
{

  // Defining member variables.

  mat = new ARSparseMatrix<ARTYPE>(np, np, pcolp, irowp, ap, nnzp);
  pA  = mat;

  this->m   = np;
  this->n   = np;

  // Checking data.

  if (check && !mat->Check()) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumNonSymMatrix::DefineMatrix");
  }

  umfpack_defaults<ARTYPE>(control);

  control[UMFPACK_PIVOT_TOLERANCE] = thresholdp;

  this->defined = true;

} // DefineMatrix (square).


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int mp, int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
             bool check, bool owner)
{

  // Defining member variables.

  mat = new ARSparseMatrix<ARTYPE>(mp, np, pcolp, irowp, ap, nnzp, '*', owner);
  pA  = mat;

  this->m  = mp;
  this->n  = np;

  // Checking data.

  if (check && !mat->Check()) {
      throw ArpackError(ArpackError::INCONSISTENT_DATA,
          "ARumNonSymMatrix::DefineMatrix");
  }

  umfpack_defaults<ARTYPE>(control);

  this->defined  = true;

} // DefineMatrix (rectangular).


template<class ARTYPE, class ARFLOAT>
inline ARumNonSymMatrix<ARTYPE, ARFLOAT>::
ARumNonSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
                 double thresholdp, bool check)
  : ARMatrix<ARTYPE>(np), mat(nullptr), AsI(nullptr)
{

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, thresholdp, check, false);

} // Long constructor (square matrix).


template<class ARTYPE, class ARFLOAT>
inline ARumNonSymMatrix<ARTYPE, ARFLOAT>::
ARumNonSymMatrix(int mp, int np, int nnzp, ARTYPE* ap,
                 int* irowp, int* pcolp) : ARMatrix<ARTYPE>(mp, np), mat(nullptr), AsI(nullptr)
{

  factored = false;
  DefineMatrix(mp, np, nnzp, ap, irowp, pcolp);

} // Long constructor (rectangular matrix).


template<class ARTYPE, class ARFLOAT>
ARumNonSymMatrix<ARTYPE, ARFLOAT>::
ARumNonSymMatrix(const std::string& name, double thresholdp, bool check)
{

  factored = false;

  ARhbMatrix<int, ARTYPE> mat;
  try {
    mat.Define(name, false);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARumNonSymMatrix");
  }

  if (mat.NCols() == mat.NRows()) {
    DefineMatrix(mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), thresholdp, check, true);
  }
  else {
    DefineMatrix(mat.NRows(), mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), check, true);
  }

} // Long constructor (Harwell-Boeing file).


template<class ARTYPE, class ARFLOAT>
ARumNonSymMatrix<ARTYPE, ARFLOAT>& ARumNonSymMatrix<ARTYPE, ARFLOAT>::
operator=(const ARumNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUNSMAT_H

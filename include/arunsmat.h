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

#include <cstddef>
#include <string>
#include "arch.h"
#include "armat.h"
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
  int     nnz;
  int*    irow;
  int*    pcol;
  int*    index;
  double  threshold;
  ARTYPE* a;
  ARTYPE* value;
  ARhbMatrix<int, ARTYPE> mat;

  bool DataOK();

  void ClearMem();

  virtual void Copy(const ARumNonSymMatrix& other);

  void SubtractAsI(ARTYPE sigma);

  void ThrowError();

 public:

  int nzeros() { return nnz; }


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

  void DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                    int* pcolp, double thresholdp = 0.1,
                    bool check = true); // Square.

  void DefineMatrix(int mp, int np, int nnzp, ARTYPE* ap,
                    int* irowp, int* pcolp);                   // Rectangular.

  ARumNonSymMatrix(): ARMatrix<ARTYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARumNonSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                   int* pcolp, double thresholdp = 0.1,
                   bool check = true);
  // Long constructor (square matrix).

  ARumNonSymMatrix(int mp, int np, int nnzp, ARTYPE* ap,
                   int* irowp, int* pcolp);
  // Long constructor (rectangular matrix).

  ARumNonSymMatrix(const std::string& name, double thresholdp = 0.1,
                   bool check = true);
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
bool ARumNonSymMatrix<ARTYPE, ARFLOAT>::DataOK()
{

  int i, j, k;

  // Checking if pcol is in ascending order.

  i = 0;
  while ((i!=this->n)&&(pcol[i]<=pcol[i+1])) i++;
  if (i!=this->n) return false;

  // Checking if irow components are in order and within bounds.

  for (i=0; i!=this->n; i++) {
    j = pcol[i];
    k = pcol[i+1]-1;
    if (j<=k) {
      if ((irow[j]<0)||(irow[k]>=this->n)) return false;
      while ((j!=k)&&(irow[j]<irow[j+1])) j++;
      if (j!=k) return false;
    }  
  }

  return true;

} // DataOK.


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::ClearMem()
{

  if (factored) {
    if (Numeric) umfpack_free_numeric<ARTYPE>(&Numeric);

    Numeric = nullptr;
  }

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
  nnz       = other.nnz;
  threshold = other.threshold;

  // Returning from here if "other" was not initialized.

  if (!this->defined) return;

  // Copying arrays with static dimension.

  for (int i = 0; i < UMFPACK_CONTROL; i++) control[i] = other.control[i];
  for (int i = 0; i < UMFPACK_INFO; i++) info[i] = other.info[i];

  // Returning from here if "other" was not factored.

  if (!factored) return;
  
} // Copy.


template<class ARTYPE, class ARFLOAT>
void ARumNonSymMatrix<ARTYPE, ARFLOAT>::SubtractAsI(ARTYPE sigma)
{

  // Subtracting sigma from diagonal elements.

  for (int i=0; i!=this->n; i++) {

    // TODO

  }

} // SubtractAsI.


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::ThrowError()
{

  int status = info[0];
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

} // ThrowError.


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

  //umfpack_symbolic(m, n, ap, ai, ax, &Symbolic, control, info);
  //umfpack_numeric(ap, ai, ax, Symbolic, &Numeric, control, info);

  //umfpack_free_symbolic<ARTYPE>(&Symbolic);

  ThrowError();

  factored = true;

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

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  void *Symbolic;

  //umfpack_symbolic(m, n, ap, ai, ax, &Symbolic, control, info);
  //umfpack_numeric(ap, ai, ax, Symbolic, &Numeric, control, info);

  //umfpack_free_symbolic<ARTYPE>(&Symbolic);

  // Handling errors.

  ThrowError();

  factored = true;

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

  // Determining w = M.v.

  for (i=0; i!=this->m; i++) w[i]=(ARTYPE)0;

  for (i=0; i!=this->n; i++) {
    t = v[i];
    for (j=pcol[i]; j!=pcol[i+1]; j++) {
      w[irow[j]] += t*a[j];
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

  for (i=0; i!=this->n; i++) {
    t = (ARTYPE)0;
    for (j=pcol[i]; j!=pcol[i+1]; j++) {
      t += v[irow[j]]*a[j];
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

  // Solving A.w = v (or AsI.w = v).

  //int status = umfpack_solve(UMFPACK_A, ap, ai, ax, w, v, Numeric, control, info);

} // MultInvv.


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
             int* pcolp, double thresholdp, bool check)
{

  // Defining member variables.

  this->m         = np;
  this->n         = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[this->n]   = nnz;
  threshold = thresholdp;

  // Checking data.

  if (check && !DataOK()) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumNonSymMatrix::DefineMatrix");
  }

  umfpack_defaults<ARTYPE>(control);

  control[UMFPACK_PIVOT_TOLERANCE] = thresholdp;

  this->defined = true;

} // DefineMatrix (square).


template<class ARTYPE, class ARFLOAT>
inline void ARumNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int mp, int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp)
{

  // Defining member variables.

  this->m        = mp;
  this->n        = np;
  nnz      = nnzp;
  a        = ap;
  irow     = irowp;
  pcol     = pcolp;
  pcol[this->n]  = nnz;

  umfpack_defaults<ARTYPE>(control);

  this->defined  = true;

} // DefineMatrix (rectangular).


template<class ARTYPE, class ARFLOAT>
inline ARumNonSymMatrix<ARTYPE, ARFLOAT>::
ARumNonSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                 int* pcolp, double thresholdp,
                 bool check): ARMatrix<ARTYPE>(np)
{

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, thresholdp, check);

} // Long constructor (square matrix).


template<class ARTYPE, class ARFLOAT>
inline ARumNonSymMatrix<ARTYPE, ARFLOAT>::
ARumNonSymMatrix(int mp, int np, int nnzp, ARTYPE* ap,
                 int* irowp, int* pcolp)             : ARMatrix<ARTYPE>(mp, np)
{

  factored = false;
  DefineMatrix(mp, np, nnzp, ap, irowp, pcolp);

} // Long constructor (rectangular matrix).


template<class ARTYPE, class ARFLOAT>
ARumNonSymMatrix<ARTYPE, ARFLOAT>::
ARumNonSymMatrix(const std::string& name, double thresholdp,
                 bool check)
{

  factored = false;

  try {
    mat.Define(name);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARumNonSymMatrix");
  }

  if (mat.NCols() == mat.NRows()) {
    DefineMatrix(mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), thresholdp, check);
  }
  else {
    DefineMatrix(mat.NRows(), mat.NCols(), mat.NonZeros(),
                 (ARTYPE*)mat.Entries(), mat.RowInd(), mat.ColPtr());
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

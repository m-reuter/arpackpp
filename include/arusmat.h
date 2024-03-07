/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUSMat.h.
   Arpack++ class ARumSymMatrix definition.

   Modified to work with Umfpack v5.??
      Martin Reuter
      Date 02/28/2013

   Arpack++ Author:
      Francisco Gomes
      
   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "aruspen.h"

#ifndef ARUSMAT_H
#define ARUSMAT_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "armat.h"
#include "arspmat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "umfpackc.h"

template<class ARTYPE> class ARumSymPencil;

template<class ARTYPE>
class ARumSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARumSymPencil<ARTYPE>;

 protected:
  
  double  control[UMFPACK_CONTROL];
  double  info[UMFPACK_INFO];
  void*   Numeric;
  bool    factored;
  char    uplo;

  // The input matrix.
  ARSparseMatrix<ARTYPE>* A;

  // In case the input matrix is triangular, UMFPACK requires the expanded matrix.
  ARSparseMatrix<ARTYPE>* Afull;

  virtual void Copy(const ARumSymMatrix& other);

  void ClearMem();

  void ExpandA();

  void SubtratcAsI(ARTYPE sigma = (ARTYPE)0);

  void Check(int status);

 private:

    // Internal matrix storing A - s I
    ARSparseMatrix<ARTYPE>* AsI;

    // Internal reference to current matrix (either A or AsI)
    ARSparseMatrix<ARTYPE>* pA;

 public:

  int nzeros() { return A->nzeros(); }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                    int* pcolp, char uplop = 'L', double thresholdp = 0.1, 
                    bool check = true, bool owner = false);

  ARumSymMatrix(): ARMatrix<ARTYPE>(), factored(false), Numeric(nullptr), A(nullptr), AsI(nullptr), Afull(nullptr)
  {
  }
  // Short constructor that does nothing.

  ARumSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                int* pcolp, char uplop = 'L', double thresholdp = 0.1,
                bool check = true);
  // Long constructor.

  ARumSymMatrix(const std::string& name, double thresholdp = 0.1,
                bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARumSymMatrix(const ARumSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumSymMatrix() { ClearMem(); }
  // Destructor.

  ARumSymMatrix& operator=(const ARumSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARumSymMatrix<ARTYPE>::ClearMem()
{

  if (factored && Numeric)
  {
    umfpack_free_numeric<ARTYPE>(&Numeric);
    Numeric = nullptr;
  }

  if (Afull && Afull != A) { delete Afull; Afull = nullptr; }

  if (A) { delete A; A = nullptr; }
  if (AsI) { delete AsI; AsI = nullptr; }

  pA = nullptr;

} // ClearMem.



template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::Copy(const ARumSymMatrix<ARTYPE>& other)
{

  // Copying very fundamental variables.
  ClearMem();

  // Copying very fundamental variables and user-defined parameters.

  this->m         = other.m;
  this->n         = other.n;
  this->defined   = other.defined;

  factored = false;

  A->Copy(*other.A);

  // Returning from here if "other" was not initialized.

  if (!this->defined) return;

  // Copying arrays with static dimension.

  for (int i = 0; i < UMFPACK_CONTROL; i++) control[i] = other.control[i];
  for (int i = 0; i < UMFPACK_INFO; i++) info[i] = other.info[i];

} // Copy.

template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::ExpandA()
{
  if (!Afull)
  {
      if (A->IsTriangular())
      {
          int ndiag = A->DiagIndices();
          int nz = 2 * (A->nzeros() - ndiag) + this->n;

          Afull = new ARSparseMatrix<ARTYPE>(this->m, this->n, nz);

          A->Expand(*Afull);
      }
      else
      {
          Afull = A;
      }
  }
}

template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::SubtratcAsI(ARTYPE sigma)
{
  if (!AsI)
  {
      ExpandA();

      int ndiag = Afull->DiagIndices();
      int nz = Afull->nzeros();

      AsI = new ARSparseMatrix<ARTYPE>(this->m, this->n, nz + this->n - ndiag);
  }

  AsI->Copy(*Afull);

  if (sigma != (ARTYPE)0)
  {
    AsI->AddDiag(-sigma);
  }

  pA = AsI;

}

template<class ARTYPE>
inline void ARumSymMatrix<ARTYPE>::Check(int status)
{

  if (status == UMFPACK_ERROR_out_of_memory)  {
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARumSymMatrix::FactorA");
  }
  else if (status == UMFPACK_ERROR_invalid_matrix)  {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumSymMatrix::FactorA");
  }
  else if (status == UMFPACK_WARNING_singular_matrix) {
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARumSymMatrix::FactorA");
  }
  else if (status != UMFPACK_OK) {
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARumSymMatrix::FactorA");
  }

} // Check.


template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::FactorA()
{

  // Quitting the function if A was not defined.
  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::FactorA");
  }

  ExpandA();

  void *Symbolic;

  auto ap = Afull->pcol();
  auto ai = Afull->irow();
  auto ax = Afull->values();

  Check(umfpack_symbolic(this->m, this->n, ap, ai, ax, &Symbolic, control, info));
  Check(umfpack_numeric(ap, ai, ax, Symbolic, &Numeric, control, info));

  umfpack_free_symbolic<ARTYPE>(&Symbolic);

  factored = true;

  pA = A;

} // FactorA.


template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.
  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::FactorAsI");
  }

  // Subtracting sigma*I from A.
  SubtratcAsI(sigma);

  // Decomposing AsI.

  void *Symbolic;

  auto ap = AsI->pcol();
  auto ai = AsI->irow();
  auto ax = AsI->values();

  Check(umfpack_symbolic(this->m, this->n, ap, ai, ax, &Symbolic, control, info));
  Check(umfpack_numeric(ap, ai, ax, Symbolic, &Numeric, control, info));

  umfpack_free_symbolic<ARTYPE>(&Symbolic);

  factored = true;

} // FactorAsI.


template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::MultMv(ARTYPE* v, ARTYPE* w)
{

  int    i,j,k;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::MultMv");
  }

  auto ax = A->values();
  auto ap = A->pcol();
  auto ai = A->irow();

  // Determining w = M.v.

  for (i = 0; i != this->m; i++) w[i] = (ARTYPE)0;

  if (uplo == 'L') {

    for (i = 0; i != this->n; i++) {
      t = v[i];
      k = ap[i];
      if (k != ap[i+1] && ai[k] == i) {
        w[i] += t*ax[k];
        k++;
      }
      for (j = k; j < ap[i+1]; j++) {
        w[ai[j]] += t*ax[j];
        w[i] += v[ai[j]]*ax[j];
      }
    }

  }
  else if (uplo == 'U') {

    for (i=  0; i != this->n; i++) {
      t = v[i];
      k = ap[i+1];
      if (k!=ap[i] && ai[k-1] == i) {
        w[i] += t*ax[k-1];
        k--;
      }
      for (j = ap[i]; j < k; j++) {
        w[ai[j]] += t*ax[j];
        w[i] += v[ai[j]]*ax[j];
      }
    }

  }
  else {

      for (i = 0; i != this->n; i++)
      {
          t = v[i];
          for (j = ap[i]; j != ap[i + 1]; j++)
          {
              w[ai[j]] += t * ax[j];
          }
      }

  }

} // MultMv.


template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARumSymMatrix::MultInvv");
  }

  auto ap = pA->pcol();
  auto ai = pA->irow();
  auto ax = pA->values();

  // Solving A.w = v (or AsI.w = v).

  int status = umfpack_solve(UMFPACK_A, ap, ai, ax, w, v, Numeric, control, info);

  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymMatrix::MultInvv");

} // MultInvv.


template<class ARTYPE>
inline void ARumSymMatrix<ARTYPE>::
DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
             int* pcolp, char uplop, double thresholdp,
             bool check, bool owner)
{
  ClearMem();

  this->m = np;
  this->n = np;

  uplo = uplop;

  A = new ARSparseMatrix<ARTYPE>(np, np, pcolp, irowp, ap, nnzp, uplop, owner);

  // Checking data.
  if (check && !A->Check()) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumSymMatrix::DefineMatrix");
  }

  umfpack_defaults<ARTYPE>(control);

  control[UMFPACK_PIVOT_TOLERANCE] = thresholdp;

  this->defined = true;

} // DefineMatrix.


template<class ARTYPE>
inline ARumSymMatrix<ARTYPE>::
ARumSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
              char uplop, double thresholdp, bool check)
    : ARMatrix<ARTYPE>(np), factored(false), Numeric(nullptr),
      A(nullptr), AsI(nullptr), Afull(nullptr)
{
  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, uplop,
               thresholdp, check);

} // Long constructor.


template<class ARTYPE>
ARumSymMatrix<ARTYPE>::
ARumSymMatrix(const std::string& file, double thresholdp, bool check)
    : ARMatrix<ARTYPE>(), factored(false), Numeric(nullptr),
      A(nullptr), AsI(nullptr), Afull(nullptr)
{
  ARhbMatrix<int, ARTYPE> mat;
  try {
    mat.Define(file, false);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARumSymMatrix");
  }

  if (mat.NCols() == mat.NRows() && mat.IsSymmetric()) {

    DefineMatrix(mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), 'L', thresholdp, check);
  }
  else {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumSymMatrix::ARluSymMatrix");
  }

} // Long constructor (Harwell-Boeing file).


template<class ARTYPE>
ARumSymMatrix<ARTYPE>& ARumSymMatrix<ARTYPE>::
operator=(const ARumSymMatrix<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSMAT_H

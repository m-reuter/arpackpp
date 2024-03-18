/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARCSMat.h.
   Arpack++ class ARchSymMatrix definition.
   (CHOLMOD wrapper)

   Author of this class:
      Martin Reuter
      Date 11/05/2012
      
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


#include "arcspen.h"

#ifndef ARCSMAT_H
#define ARCSMAT_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "armat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "cholmodc.h"

template<class ARTYPE> class ARchSymPencil;

template<class ARTYPE>
class ARchSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARchSymPencil<ARTYPE>;

 protected:

  bool    factored;
  char    uplo;
  int     nnz;
  int*    irow;
  int*    pcol;
  ARTYPE* a;
  cholmod_common c ;
  cholmod_sparse *A ; 
  cholmod_factor *L ; 
   
  bool DataOK();

  virtual void Copy(const ARchSymMatrix& other);

  void ClearMem();

 public:

  int nzeros() { return nnz; }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                    int* pcolp, char uplop = 'L', bool check = true);

  ARchSymMatrix(): ARMatrix<ARTYPE>(), factored(false) { cholmod_start(&c); }
  // Short constructor that does nothing.

  ARchSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                int* pcolp, char uplop = 'L', bool check = true);
  // Long constructor.

  ARchSymMatrix(const std::string& name, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARchSymMatrix(const ARchSymMatrix& other) { cholmod_start(&c); Copy(other); }
  // Copy constructor.

  virtual ~ARchSymMatrix() { ClearMem(); cholmod_finish(&c); }
  // Destructor.

  ARchSymMatrix& operator=(const ARchSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARchSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
bool ARchSymMatrix<ARTYPE>::DataOK()
{

  return cholmod_check_sparse(A, &c) == 1 /* TRUE */;

} // DataOK.


template<class ARTYPE>
void ARchSymMatrix<ARTYPE>::ClearMem()
{

  if (factored) {
    cholmod_free_factor (&L, &c) ;
  }
  if (this->defined) {
    //cholmod_free_sparse (&A, &c);

    free(A); // don't delete data in A as it came from external
    A = nullptr;
  }

} // ClearMem.



template<class ARTYPE>
inline void ARchSymMatrix<ARTYPE>::Copy(const ARchSymMatrix<ARTYPE>& other)
{

  // Copying very fundamental variables.
  ClearMem();

  this->defined   = other.defined;
  // Returning from here if "other" was not initialized.
  if (!this->defined) return;

  this->n = other.n;
  factored  = other.factored;
  uplo = other.uplo;
  nnz  = other.nnz;
  irow = other.irow;
  pcol = other.pcol;
  a = other.a;
  //c = other.c;
   
  A = cholmod_copy_sparse(other.A,&c);

  if (L) cholmod_free_factor(&L,&c);
  if (factored)
    L = cholmod_copy_factor(other.L,&c);

} // Copy.



template<class ARTYPE>
void ARchSymMatrix<ARTYPE>::FactorA()
{
  int info;

  // Quitting the function if A was not defined.
  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARchSymMatrix::FactorA");
  }

  // Deleting previous versions of L.
  if (factored) {
    cholmod_free_factor (&L, &c) ;
  }
  
  L = cholmod_analyze (A, &c) ;
  info = cholmod_factorize (A, L, &c) ;  
  

  factored = (info != 0);
  
  if (c.status != CHOLMOD_OK)
  {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARchSymMatrix::FactorA");
    
    factored = false;
  }
  
// 
//   // Handling errors.
// 
//   if (info < 0)  {        // Illegal argument.
//     throw ArpackError(ArpackError::PARAMETER_ERROR,
//                       "ARchSymMatrix::FactorA");
//   }
//   else if (info > this->n) {    // Memory is not sufficient.
//     throw ArpackError(ArpackError::MEMORY_OVERFLOW,
//                       "ARchSymMatrix::FactorA");
//   }
//   else if (info > 0) {   // Matrix is singular.
//     throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
//                       "ARchSymMatrix::FactorA");
//   }

} // FactorA.


template<class ARTYPE>
void ARchSymMatrix<ARTYPE>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.
  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARchSymMatrix::FactorAsI");
  }

  // Deleting previous versions of L.
  if (factored) {
    cholmod_free_factor (&L, &c);
  }  

  // Factorizing A-sigma*I
  double sigma2[2] = { -sigma, 0.0 };
  L = cholmod_analyze (A, &c);
  int info = cholmod_factorize_p (A, sigma2, nullptr, 0, L, &c);

  factored = (info != 0);
  
  if (c.status != CHOLMOD_OK)
  {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARchSymMatrix::FactorAsI");
    
    factored = false;
  }

} // FactorAsI.


template<class ARTYPE>
void ARchSymMatrix<ARTYPE>::MultMv(ARTYPE* v, ARTYPE* w)
{
  int    i, j, k;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARchSymMatrix::MultMv");
  }

  // Determining w = M.v.

  for (i=0; i!=this->m; i++) w[i]=(ARTYPE)0;

  if (uplo == 'U') {

    for (i=0; i!=this->n; i++) {
      t = v[i];
      k = pcol[i+1];
      if ((k!=pcol[i])&&(irow[k-1]==i)) {
        w[i] += t*a[k-1];
        k--;
      }
      for (j=pcol[i]; j<k; j++) {
        w[irow[j]] += t*a[j];
        w[i] += v[irow[j]]*a[j];
      }
    }

  }
  else {

    for (i=0; i!=this->n; i++) {
      t = v[i];
      k = pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) {
        w[i] += t*a[k];
        k++;
      }
      for (j=k; j<pcol[i+1]; j++) {
        w[irow[j]] += t*a[j];
        w[i] += v[irow[j]]*a[j];
      }
    }

  }

} // MultMv.


template<class ARTYPE>
void ARchSymMatrix<ARTYPE>::MultInvv(ARTYPE* v, ARTYPE* w)
{
  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARchSymMatrix::MultInvv");
  }

  // Solving A.w = v (or AsI.w = v).
  
  //create b from v (data is not copied!!)
  cholmod_dense *b = CholmodCreateDense(this->n, 1, v);

  cholmod_dense *x = cholmod_solve (CHOLMOD_A, L, b, &c) ;

  CholmodGetDenseData(x, this->n, w);

  free(b);
  cholmod_free_dense(&x, &c);

} // MultInvv.


template<class ARTYPE>
inline void ARchSymMatrix<ARTYPE>::
DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
             char uplop, bool check)
{

  this->m   = np;
  this->n   = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[this->n]   = nnz;
  uplo      = uplop;

  // Creating cholmod_sparse A.
  A = CholmodCreateSparse(this->n, this->n, nnz, a, irow, pcol, uplo);

  this->defined = true;

  // Checking data.
  if (check && !DataOK()) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARchSymMatrix::DefineMatrix");
  }

} // DefineMatrix.


template<class ARTYPE>
inline ARchSymMatrix<ARTYPE>::
ARchSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
              int* pcolp, char uplop, bool check) : ARMatrix<ARTYPE>(np)
{
 cholmod_start (&c) ;

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, uplop, check);

} // Long constructor.


template<class ARTYPE>
ARchSymMatrix<ARTYPE>::
ARchSymMatrix(const std::string& file, bool check)
{
 cholmod_start (&c) ;

  factored = false;

  ARhbMatrix<int, ARTYPE> mat;

  try {
    mat.Define(file, false);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARchSymMatrix");
  }

  if ((mat.NCols() == mat.NRows()) && (mat.IsSymmetric())) {

    DefineMatrix(mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), 'L', check);
  }
  else {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARchSymMatrix::ARchSymMatrix");
  }
} // Long constructor (Harwell-Boeing file).


template<class ARTYPE>
ARchSymMatrix<ARTYPE>& ARchSymMatrix<ARTYPE>::
operator=(const ARchSymMatrix<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARCSMAT_H

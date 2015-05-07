/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE AREig.h.
   Example functions AREig. This function returns eigenvalues, EigVal, and
   eigenvectors, EigVec, of real symmetric, real nonsymmetric and complex
   problems in various modes using ARluSymStdEig, ARluNonSymStdEig and
   ARluNonSymGenEig classes.
   The SuperLU package is employed to solve the linear systems that appear
   when the shift-and-invert spectral transformation is being used.

   There are eighteen different versions of AREig, as shown below. The type
   and the meaning of each AREig parameter is briefly described in section
   II. For a complete description of all parameters, see the Appendix of
   "ARPACK++: a c++ implementation of ARPACK eigenvalue package."

   I) How to call AREig:

      AREig is a c++ overloaded function, which means that there are
      several definitions for this single function. Each definition is
      related to a different problem and a different kind of data. All
      twenty six available AREig functions are listed below. They are
      divided according to the problem type (standard or generalized),
      the type of the involved matrices (real or complex), the
      computational mode being used (regular or shift-and-invert) and
      the desired output data (only eigenvalues or eigenvalues and
      eigenvectors).

      1) Real symmetric standard problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, uplo, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, uplo,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

         b) Real shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, uplo, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, uplo, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

      2) Real symmetric generalized problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, uplo, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, uplo, nev, which,
                          ncv, tol, maxit, resid, AutoShift)

         b) Shift-and-invert, buckling and Cayley modes:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, uplo, InvertMode, sigma, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA, nnzB,
                          B, irowB, pcolB, uplo, InvertMode, sigma, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

      3) Real nonsymmetric standard problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnz, A, irow, pcol, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

         b) Real shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnz, A, irow, pcol, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol,
                          sigma, nev, which, ncv, tol, maxit, resid, AutoShift)

      4) Real nonsymmetric generalized problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnzA, A, irowA,
                          pcolA, nnzB, B, irowB, pcolB, nev, which, ncv,
                          tol, maxit, resid, AutoShift)

         b) Real shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, sigma, nev, which, ncv,
                          tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnzA, A, irowA,
                          pcolA, nnzB, B, irowB, pcolB, sigma, nev, which,
                          ncv, tol, maxit, resid, AutoShift)

         c) Complex shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, part, sigmaR, SigmaI,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, part, sigmaR, SigmaI,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

      5) Complex standard problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

         b) Shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

      6) Complex generalized problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, nev, which, ncv, tol, maxit,
                          resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

         b) Shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, sigma, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, sigma, nev, which,
                          ncv, tol, maxit, resid, AutoShift)


   II) AREig parameters:

   Each AREig function contains a variable number of parameters, because
   some parameters are not mandatory. Optional parameters should only be
   defined by the used if the default values supplied by the function are
   not suitable.
   In the description of AREig parameters, the following notation is used
   to describe eigenvalue problems:
   a) standard problems:    A*EigVec = Eigvec*EigVal,
   b) generalized problems: A*EigVec = B*EigVec*EigVal.

      1) Compulsory input parameters:

      int n             Dimension of the problem.

      int nnz           Number of nonzero elements in matrix A.

      int nnzA          Same as nnz.

      int nnzB          Number of nonzero elements in matrix B.

      TYPE A[]          Array of nonzero elements in matrix A.
                        TYPE must be one of "float", "double",
                        "arcomplex<float>" or "arcomplex<double>".

      TYPE B[]          Array of nonzero elements in matrix B.
                        TYPE must be one of "float", "double",
                        "arcomplex<float>" or "arcomplex<double>".

      int irow[]        Array of row indices of the nonzero elements in A.

      int irowA[]       Same as irow.

      int irowB[]       Array of row indices of the nonzero elements in B.

      int pcol[]        Array of pointers to the beginning of columns in
                        A and irow. pcol must have n+1 elements and the
                        last element must be nnz.

      int pcolA[]       Same as pcol.

      int pcolB[]       Array of pointers to the beginning of columns in
                        B and irowB. pcol must have n+1 elements and the
                        last element must be nnzB.

      char uplo         A parameter used only if the problem is symmetric.
                        uplo indicates whether the lower (uplo = "L") or the
                        upper triangular (uplo = "U") part of A (and also
                        B, if the problem is a generalized one) is being
                        supplied by the user.

      char InvertMode   Spectral transformation used when solving symmetric
                        generalized problems. To use the shift and invert
                        mode, the user must set InvertMode to 'S'.  Buckling
                        and Cayley modes are represented by 'B' and 'C',
                        respectively.

      TYPE sigma        The shift (when shift-and invert mode is being used).
                        TYPE must be one of "float" or "double" if the
                        problem is nonsymmetric and sigma is real, or one of
                        "arcomplex<float>" or "arcomplex<double>" if the
                        problem is complex.

      ARFLOAT sigmaR    Real part of the shift if the problem is real but
                        the shift is complex. ARFLOAT must be one of "float"
                        or "double".

      ARFLOAT sigmaI    Imaginary part of the shift if the problem is real
                        but the shift is complex. ARFLOAT must be one of
                        "float" or "double".

      char part         A parameter that characterizes which part (real or
                        imaginary) of the vector y = OP*x will be used by
                        ARPACK++ when the problem is real but the shift is
                        complex. "part" must be set to one of 'R' (real
                        part) or 'I' (imaginary part).

      int  nev          Number of eigenvalues to be computed.


      2) Optional input parameters:

      const std::string& which
                        A parameter thar specifies which of the Ritz values
                        are to be computed. "which" must be set to one of:
                        LM: to find eigenvalues with largest magnitude;
                        SM: to find eigenvalues with smallest magnitude;
                        LR: to find eigenvalues with largest real part;
                        SR: to find eigenvalues with smallest real part;
                        LI: to find eigenvalues with largest imaginary part;
                        SI: to find eigenvalues with smallest imaginary part.
                        Default: LM.

      int  ncv          Number of Arnoldi vectors generated at each
                        iteration of ARPACK. Default: 2*nev+1.

      ARFLOAT tol       Stopping criterion (relative accuracy of Ritz
                        values). ARFLOAT must be one of "float" or "double".
                        Default: machine precision.

      int  maxit        Maximum number of Arnoldi update iterations allowed.
                        Default: 100*nev.

      TYPE* resid       A pointer to an array that contains the initail
                        vector. Default: a random vector.
                        TYPE must be one of "float", "double",
                        "arcomplex<float>" or "arcomplex<double>".

      bool AutoShift    A parameter that indicates if exact shifts for
                        implicit restarting of the Arnoldi method are to
                        be generated internally by ARPACK++ or shifts are
                        being supplied by the user. Default: true (exact
                        shifts are being used).

      3) Output parameters:

      TYPE  EigVal[]    A vector that contains the "converged" eigenvalues
                        when the problem is symmetric or complex. TYPE
                        must be one of float, double, "arcomplex<float>" or
                        "arcomplex<double>".

      ARFLOAT EigValR[] A vector that contains the real part of the
                        "converged" eigenvalues (when the problem is real).
                        ARFLOAT must be one of "float" or "double".

      ARFLOAT EigValI[] A vector that contains the imaginary part of the
                        "converged" eigenvalues (when the problem is real).
                        ARFLOAT must be one of "float" or "double".

      TYPE  EigVec[]    A vector that stores all "converged" eigenvectors
                        consecutively. For real problems, complex eigenvectors
                        are given as two consecutive vectors. The first
                        contains the real part of the eigenvector, while
                        the imaginary part is stored in the second vector.
                        TYPE must be one of "float", "double",
                        "arcomplex<float>" or "arcomplex<double>".


   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef AREIG_H
#define AREIG_H

#include <string>
#include "arcomp.h"
#include "arlsmat.h"
#include "arlnsmat.h"
#include "arlssym.h"
#include "arlgsym.h"
#include "arlsnsym.h"
#include "arlgnsym.h"
#include "arlscomp.h"
#include "arlgcomp.h"


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], int n, int nnz, arcomplex<ARFLOAT> A[],
          int irow[], int pcol[], int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, arcomplex<ARFLOAT>* resid = NULL,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<ARFLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex standard problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], arcomplex<ARFLOAT> EigVec[], int n,
          int nnz, arcomplex<ARFLOAT> A[], int irow[], int pcol[],
          int nev, const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, arcomplex<ARFLOAT>* resid = NULL, 
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<ARFLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex standard problem, values and vectors, regular mode.


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], int n, int nnz, arcomplex<ARFLOAT> A[],
          int irow[], int pcol[], arcomplex<ARFLOAT> sigma, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0,
          arcomplex<ARFLOAT>* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<ARFLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex standard problem, only eigenvalues, shift-and-invert.


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], arcomplex<ARFLOAT> EigVec[], int n,
          int nnz, arcomplex<ARFLOAT> A[], int irow[], int pcol[],
          arcomplex<ARFLOAT> sigma, int nev, const std::string& which = "LM",
          int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0,
          arcomplex<ARFLOAT>* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<ARFLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex standard problem, values and vectors, shift-and-invert.


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], int n, int nnzA,
          arcomplex<ARFLOAT> A[], int irowA[], int pcolA[], int nnzB,
          arcomplex<ARFLOAT> B[], int irowB[], int pcolB[], int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, arcomplex<ARFLOAT>* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixA(n,nnzA,A,irowA,pcolA);
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixB(n,nnzB,B,irowB,pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<ARFLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex generalized problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], arcomplex<ARFLOAT> EigVec[], int n,
          int nnzA, arcomplex<ARFLOAT> A[], int irowA[], int pcolA[],
          int nnzB, arcomplex<ARFLOAT> B[], int irowB[], int pcolB[],
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, arcomplex<ARFLOAT>* resid = NULL, 
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixA(n,nnzA,A,irowA,pcolA);
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixB(n,nnzB,B,irowB,pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<ARFLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex generalized problem, values and vectors, regular mode.


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], int n, int nnzA, arcomplex<ARFLOAT> A[],
          int irowA[], int pcolA[], int nnzB, arcomplex<ARFLOAT> B[],
          int irowB[], int pcolB[], arcomplex<ARFLOAT> sigma, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, arcomplex<ARFLOAT>* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixA(n,nnzA,A,irowA,pcolA);
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixB(n,nnzB,B,irowB,pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<ARFLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex generalized problem, only eigenvalues, shift-and-invert mode.


template <class ARFLOAT>
int AREig(arcomplex<ARFLOAT> EigVal[], arcomplex<ARFLOAT> EigVec[], int n,
          int nnzA, arcomplex<ARFLOAT> A[], int irowA[], int pcolA[],
          int nnzB, arcomplex<ARFLOAT> B[], int irowB[], int pcolB[],
          arcomplex<ARFLOAT> sigma, int nev, const std::string& which = "LM",
          int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0,
          arcomplex<ARFLOAT>* resid = NULL, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixA(n,nnzA,A,irowA,pcolA);
  ARluNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> matrixB(n,nnzB,B,irowB,pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<ARFLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex generalized problem, values and vectors, shift-and-invert mode.


template <class ARFLOAT>
int AREig(double EigValR[], ARFLOAT EigValI[], int n, int nnz,
          ARFLOAT A[], int irow[], int pcol[], int nev,
          const std::string& which = "LM",
          int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0, 
          ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<ARFLOAT> prob(nev, matrix, which, ncv, tol,
                                 maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(float EigValR[], ARFLOAT EigValI[], int n, int nnz,
          ARFLOAT A[], int irow[], int pcol[], int nev,
          const std::string& which = "LM",
          int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0, 
          ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<ARFLOAT> prob(nev, matrix, which, ncv, tol,
                                 maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigValR[], ARFLOAT EigValI[], ARFLOAT EigVec[], int n, 
          int nnz, ARFLOAT A[], int irow[], int pcol[], int nev, 
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0, 
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<ARFLOAT> prob(nev, matrix, which, ncv, tol,
                                 maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric standard problem, values and vectors, regular mode.


template <class ARFLOAT>
int AREig(double EigValR[], ARFLOAT EigValI[], int n, int nnz,
          ARFLOAT A[], int irow[], int pcol[], ARFLOAT sigma, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<ARFLOAT> prob(nev, matrix, sigma, which, ncv,
                                 tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, shift-and-invert.


template <class ARFLOAT>
int AREig(float EigValR[], ARFLOAT EigValI[], int n, int nnz,
          ARFLOAT A[], int irow[], int pcol[], ARFLOAT sigma, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<ARFLOAT> prob(nev, matrix, sigma, which, ncv,
                                 tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, shift-and-invert.


template <class ARFLOAT>
int AREig(ARFLOAT EigValR[], ARFLOAT EigValI[], ARFLOAT EigVec[], int n, 
          int nnz, ARFLOAT A[], int irow[], int pcol[], ARFLOAT sigma, 
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<ARFLOAT> prob(nev, matrix, sigma, which, ncv,
                                 tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric standard problem, values and vectors, shift-and-invert.


template <class ARFLOAT>
int AREig(double EigValR[], ARFLOAT EigValI[], int n, int nnzA,
          ARFLOAT A[], int irowA[], int pcolA[], int nnzB,
          ARFLOAT B[], int irowB[], int pcolB[], int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, which,
                                 ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(float EigValR[], ARFLOAT EigValI[], int n, int nnzA,
          ARFLOAT A[], int irowA[], int pcolA[], int nnzB,
          ARFLOAT B[], int irowB[], int pcolB[], int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, which,
                                 ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigValR[], ARFLOAT EigValI[], ARFLOAT EigVec[], int n,
          int nnzA, ARFLOAT A[], int irowA[], int pcolA[],
          int nnzB, ARFLOAT B[], int irowB[], int pcolB[],
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, which,
                                 ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors, regular mode.


template <class ARFLOAT>
int AREig(double EigValR[], ARFLOAT EigValI[], int n, int nnzA,
          ARFLOAT A[], int irowA[], int pcolA[], int nnzB,
          ARFLOAT B[], int irowB[], int pcolB[], ARFLOAT sigma,
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, sigma, which,
                                 ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // real shift-and-invert mode.


template <class ARFLOAT>
int AREig(float EigValR[], ARFLOAT EigValI[], int n, int nnzA,
          ARFLOAT A[], int irowA[], int pcolA[], int nnzB,
          ARFLOAT B[], int irowB[], int pcolB[], ARFLOAT sigma,
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, sigma, which,
                                 ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // real shift-and-invert mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigValR[], ARFLOAT EigValI[], ARFLOAT EigVec[], int n,
          int nnzA, ARFLOAT A[], int irowA[], int pcolA[], int nnzB,
          ARFLOAT B[], int irowB[], int pcolB[], ARFLOAT sigma, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, sigma, which,
                                 ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors,
  // real shift-and-invert mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigValR[], ARFLOAT EigValI[], int n, int nnzA, ARFLOAT A[],
          int irowA[], int pcolA[], int nnzB, ARFLOAT B[], int irowB[],
          int pcolB[], char part, ARFLOAT sigmaR, ARFLOAT sigmaI,
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, part,
                                 sigmaR, sigmaI, which, ncv,
                                 tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // complex shift-and-invert mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigValR[], ARFLOAT EigValI[], ARFLOAT EigVec[], int n, 
          int nnzA, ARFLOAT A[], int irowA[], int pcolA[], int nnzB, 
          ARFLOAT B[], int irowB[], int pcolB[], char part, ARFLOAT sigmaR, 
          ARFLOAT sigmaI, int nev, const std::string& which = "LM",
          int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0,
          ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<ARFLOAT, ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, part,
                                 sigmaR, sigmaI, which, ncv,
                                 tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors,
  // complex shift-and-invert mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], int n, int nnz, ARFLOAT A[], int irow[],
          int pcol[], char uplo, int nev, const std::string& which = "LM",
          int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0,
          ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<ARFLOAT> prob(nev, matrix, which, ncv, tol,
                              maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric standard problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], ARFLOAT EigVec[], int n, int nnz, ARFLOAT A[],
          int irow[], int pcol[], char uplo, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<ARFLOAT> prob(nev, matrix, which, ncv, tol,
                              maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric standard problem, values and vectors, regular mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], int n, int nnz, ARFLOAT A[], int irow[],
          int pcol[], char uplo, ARFLOAT sigma, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<ARFLOAT> prob(nev, matrix, sigma, which, ncv,
                              tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric standard problem, only eigenvalues, shift-and-invert.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], ARFLOAT EigVec[], int n, int nnz, ARFLOAT A[],
          int irow[], int pcol[], char uplo, ARFLOAT sigma,
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<ARFLOAT> prob(nev, matrix, sigma, which, ncv,
                              tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric standard problem, values and vectors, shift-and-invert.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], int n, int nnzA, ARFLOAT A[], int irowA[],
          int pcolA[], int nnzB, ARFLOAT B[], int irowB[], int pcolB[],
          char uplo, int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, which,
                              ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric generalized problem, only eigenvalues, regular mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], ARFLOAT EigVec[], int n, int nnzA, ARFLOAT A[],
          int irowA[], int pcolA[], int nnzB, ARFLOAT B[], int irowB[],
          int pcolB[], char uplo, int nev, const std::string& which = "LM",
          int ncv = 0, ARFLOAT tol = 0.0, int maxit = 0,
          ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<ARFLOAT> prob(nev, matrixA, matrixB, which,
                              ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors, regular mode.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], int n, int nnzA, ARFLOAT A[], int irowA[],
          int pcolA[], int nnzB, ARFLOAT B[], int irowB[], int pcolB[],
          char uplo, char InvertMode, ARFLOAT sigma, int nev,
          const std::string& which = "LM", int ncv = 0, ARFLOAT tol = 0.0,
          int maxit = 0, ARFLOAT* resid = NULL, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<ARFLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                              which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric generalized problem, only eigenvalues,
  // shift-and-invert, buckling and Cayley modes.


template <class ARFLOAT>
int AREig(ARFLOAT EigVal[], ARFLOAT EigVec[], int n, int nnzA, ARFLOAT A[],
          int irowA[], int pcolA[], int nnzB, ARFLOAT B[], int irowB[],
          int pcolB[], char uplo, char InvertMode, ARFLOAT sigma,
          int nev, const std::string& which = "LM", int ncv = 0,
          ARFLOAT tol = 0.0, int maxit = 0, ARFLOAT* resid = NULL,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<ARFLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<ARFLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<ARFLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                              which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors,
  // shift-and-invert, buckling and Cayley modes.


#endif // AREIG_H

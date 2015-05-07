/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LSVD.cc.
   Example program that illustrates how to determine the truncated SVD
   of a matrix using the ARSymStdEig class.

   1) Problem description:

      In this example, Arpack++ is called to solve the symmetric problem:

                              | 0  A |*y = sigma*y,
                              | A' 0 |

      where A is an m by n real matrix.
      This problem can be used to obtain the decomposition A = U*S*V'.
      The positive eigenvalues of this problem are the singular values 
      of A (the eigenvalues come in pairs, the negative eigenvalues have
      the same magnitude of the positive ones and can be discarded). 
      The columns of U can be extracted from the first m components of
      the eigenvectors y, while the columns of V can be
      extracted from the the remaining n components.

   2) Data structure used to represent the matrix:

      {nnzA, irowA, pcolA, valA}: matrix A data in CSC format.

   3) Included header files:

      File             Contents
      -----------      --------------------------------------------
      lnmatrxv.h       RectangularMatrix, a function that generates
                       matrix A in CSC format.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arssym.h         The ARSymStdEig class definition.
      lsvdsol.h        The Solution function definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arssym.h"
#include "lnmatrxv.h"
#include "arlnsmat.h"
#include "lsvdsol.h"


int main()
{

  // Defining variables;

  int     m;          // Number of rows in A.
  int     n;          // Number of columns in A.
  int     nnz;        // Number of nonzero elements in A.
  int*    irow;       // pointer to an array that stores the row
                      // indices of the nonzeros in A.
  int*    pcol;       // pointer to an array of pointers to the
                      // beginning of each column of A in valA.
  double* valA;       // pointer to an array that stores the
                      // nonzero elements of A.

  // Creating a rectangular matrix with m = 200 and n = 100.

  n = 100;
  RectangularMatrix(n, m, nnz, valA, irow, pcol);

  // Using ARluNonSymMatrix to store matrix information and to
  // perform the product OP*x (no LU decomposition is performed).

  ARluNonSymMatrix<double, double> A(m, n, nnz, valA, irow, pcol);

  // Defining what we need: the four eigenvalues with largest
  // algebraic value.

  ARSymStdEig<double, ARluNonSymMatrix<double, double> >
    dprob(m+n, 5L, &A, &ARluNonSymMatrix<double, double>::Mult0MMt0v, 
          "LA", 20L);

  // Finding eigenvalues.

  dprob.FindEigenvectors();

  // Printing the solution.

  Solution(A, dprob);

} // main.


/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE HSymStd.cc.
   Example program that illustrates how to solve a real
   symmetric standard eigenvalue problem using the
   ARluSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in
      regular or shift and invert mode, where A is read
      from a file in Harwell-Boeing format.

   2) Included header files:

      File             Contents
      -----------      -------------------------------------------
      arerror.h        The ArpackError class definition.
      arlsmat.h        The ARluSymMatrix class definition.
      arlssym.h        The ARluSymStdEig class definition.
      lsymsol.h        The Solution function.

   3) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include <iostream>
#include <string>
#include "arerror.h"
#include "arlsmat.h"
#include "arlssym.h"
#include "lsymsol.h"


void PrintHelp()
/*
  Prints information about hsymstd usage.
*/
{

  std::cout << "ARPACK++ version 1.2 feb 2000" << std::endl;
  std::cout << "hsymstd: a standard symmetric eigenvalue problems solver" << std::endl;
  std::cout << "usage:   hsymstd [parameters] file" << std::endl;
  std::cout << "parameters:" << std::endl;
  std::cout << "      -n (number of desired eigenvalues)" << std::endl;
  std::cout << "      -c (number of Arnoldi vectors per iteration)" << std::endl;
  std::cout << "      -l (maximum number of iterations)" << std::endl;
  std::cout << "      -s (shift)" << std::endl;
  std::cout << "      -t (stopping criterion)" << std::endl;
  std::cout << "      -u (LU pivot threshold)" << std::endl;
  std::cout << "      -o (column ordering for factorization)" << std::endl;
  std::cout << "      -w (desired portion of the spectrum. " << std::endl;
  std::cout << "          acceptable values: LM, SM, LA, SA, BE)" << std::endl;
  std::cout << std::endl;

} // PrintHelp.


bool ReadParameters(int n, char* v[], int &nev, int &ncv, int &maxit,
                    int &order, bool &shift, double &sigma, double &tol,
                    double &thresh, std::string &which, std::string &file)
/*
  Reads parameters from the command line.
*/
{

  int  i;
  bool ok;

  // Defining default values.  

  nev     = 5;
  ncv     = 0;
  maxit   = 0;
  order   = 2;
  shift   = false;
  sigma   = 0.0;
  tol     = 0.0;
  thresh  = 0.1;
  which   = "LM";
  file    = " ";
  ok      = true;

  // Returning if the number of parameters is even.

  if (n%2) {
    ok = false;
  }
  else {

    // Reading parameters.

    i = 1;
    while ((i<(n-1)) && (ok)) {
    
      if ((v[i][0] != '-') || (strlen(v[i]) != 2)) {
        ok = false;
        i += 2;
      }
      else {
        switch (v[i++][1]) {
        case 'n':
          nev = atoi(v[i++]);
          break;
        case 's':
          sigma = atof(v[i++]);
          shift = true;
          break;
        case 'w':
          which = v[i++];
          break;
        case 'c':
          ncv = atoi(v[i++]);
          break;
        case 't':
          tol = atof(v[i++]);
          break;
        case 'l':
          maxit = atoi(v[i++]);
          break;
        case 'o':
          order = atoi(v[i++]);
          break;
        case 'u':
          thresh = atof(v[i++]);
          break;
        default :
          std::cout << "unrecognized parameter: -" << v[i-1][1] << std::endl;
          ok = false; 
        }
      }
    }
  }

  if (ok) {
    file = v[i];
  }
  else {
    PrintHelp();
  }

  return ok;

} // ReadParameters.


int main(int argc, char* argv[])
{

  // Defining variables.

  int    nev;
  int    ncv;
  int    maxit;
  int    order;
  bool   shift;
  double sigma;
  double tol;
  double thresh;
  std::string  which;
  std::string  file;

  // Reading parameters.

  if (!ReadParameters(argc, argv, nev, ncv, maxit, order,
                      shift, sigma, tol, thresh, which, file)) {
    return 0;
  }

  // Reading and storing the matrix.

  ARluSymMatrix<double> matrix(file, thresh, order);

  // Defining the eigenvalue problem.

  ARluSymStdEig<double> dprob(nev, matrix, which, ncv, tol, maxit);

  // Defining the shift.

  if (shift) dprob.SetShiftInvertMode(sigma);

  // Finding eigenvalues and eigenvectors.

  try {
    dprob.FindEigenvectors();
  }
  catch (ArpackError) { return 0; }

  // Printing solution.

  Solution(matrix, dprob);

} // main.


# arpackpp (ARPACK++)
# 
[![Build Status](https://travis-ci.org/m-reuter/arpackpp.svg?branch=master)](https://travis-ci.org/m-reuter/arpackpp)


## Introduction
## 
Arpackpp is a C++ interface to the ARPACK Fortran package, which
implements the implicit restarted Arnoldi method for iteratively solving
large-scale sparse eigenvalue problems.

Arpackpp is a collection of classes (C++ headers and examples) that
offers C++ programmers an interface to ARPACK. Furthermore, it interfaces
with LAPACK, SuperLU, Cholmod and UMFPACK to incorporate efficient
matrix solvers. Arpackpp preserves the full capability, performance,
accuracy and low memory requirements of the ARPACK Fortan package, but
takes advantage of the C++ object-oriented programming environment.

This GitHub project is designed to provide a common maintained version
of arpackpp. It is derived from the orignial package (ARPACK++ Version
1.2. by Gomes and Sorensen), which has not been actively maintained for
many years. Several updates have been included in this version (some of
them were previously hosted as patches at
http://reuter.mit.edu/software/arpackpatch/ ). This GitHub repository is
designed to collect fixes and updates (e.g. to more recent or future
releases of the involved libraries). Please consider contributing (see
todo list below).


## Features:
## 
Features of original ARPACK++ package:

- Friendly interface that hides the complicated reverese communication
  interface of the Fortran Arpack package from the user.
  
- Easy interface using  matrices and vectors via the Standard Template
  Library (STL). 

- Provides an interface between ARPACK and solvers in SuperLU, LAPACK,
  UMFPACK, and CHOLMOD to solve eigenvalue problems (specifically shift
  invert methods). 

- Use of templates for optimal performance.


Additional features of this GitHub arpackpp package:

- CMake support for building the examples

- Examples build on Mac OSX using CLang 

- Install scripts for getting and building the dependencies of examples 

- Support for SuperLU 5.0 http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ 

- Added initial support for CHOLMOD (for symmetric real problems) 

- Updated UMFPACK sym integration with SuiteSparse 
  http://faculty.cse.tamu.edu/davis/suitesparse.html

- Fixed ARPACK++1.2 to run with g++ 4.4.6 and SuperLU 4.3 (patch see
  here: http://reuter.mit.edu/software/arpackpatch/ )


## TODO
## 
- CMake: get rid of globbing and specify individual files, also add some testing
- UMFPACK complex examples do not build (need update like sym) 
- CHOLMOD complex examples not included (implement similar to real sym)
- Update documentation (install) to cover more scenarios (APT, Homebrew)


## Files
## 
1. Files included in the main directory:

  1. README.md:

      This file.

  2. INSTALL.md:

      Compile and install notes.

  3. Makefile.inc (historic):

      An include file used to compile arpackpp examples. You must change
      some directories and machine-dependent directives contained into
      this file prior to compiling the examples. See the description of
      the "makefiles" directory below.

  4. CmakeLists.txt:

      A Cmake file to compile arpackpp examples.

  5. install-*.sh

      Shell scripts to download and install dependencies into a local
      ./external directory. Some dependencies can also be installed via
      a package-manager on your system.


2. arpackpp subdirectories:

  1. makefiles (historic)

      This directory contains example Makefile.inc include files for
      some platforms. Choose one and copy it onto the
      arpackpp/Makefile.inc file.

  2. include:

      The directory that contains arpackpp library, i.e., all header
      files that define arpackpp class templates.

  3. examples:

      The directory where all arpackpp examples can be found. These
      examples are intended to illustrate how to use and compile
      arpackpp classes and are divided according to the type of problem
      being solved and also the kind of information that the user is
      supposed to supply. Look at the examples/README file for further
      information.

      Note: additional header files are contained in examples/matrices
      and examples/matprod that are needed to build examples or your own
      code!

  4. doc:

      The directory that contains a the arpackpp user's manual and some
      instructions on how to install the libraries required by arpackpp.

## Dependencies
## 
- LAPACK

- BLAS (e.g. OpenBLAS https://github.com/xianyi/OpenBLAS.git )

- ARPACK (arpack-ng https://github.com/opencollab/arpack-ng.git )

For specific operations only, any of these: 

- SuperLU 
  http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.0.tar.gz

- UMFPACK 
  http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.5.tar.gz 

- CHOLMOD 
  http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.5.tar.gz 


1. ARPACK (fortran):

   Arpackpp is a C++ interface to ARPACK fortran code, so the original
   ARPACK library must be installed prior to using the C++ version. A
   mainted package (arpack new generation) can be obtained via this
   GitHub repository (see also install-arpack-ng.sh):

   https://github.com/opencollab/arpack-ng

2. BLAS and LAPACK (fortran versions):

   BLAS and LAPACK routines required by ARPACK fortran code are
   distributed along with the software. However, some arpackpp examples
   require routines from these libraries that are not included in the
   ARPACK distribution, so it is recommended to install BLAS and LAPACK
   before compiling the examples. Besides, you should use
   vendor-optimized versions of these libraries if they are available.
   E.g. OpenBLAS is available via this GitHub repository (see also
   install-openblas.sh):

   https://github.com/xianyi/OpenBLAS

3. SUPERLU:

   Some ARPACK++ classes call SUPERLU library functions to solve
   eigenvalue problems that require complex or real (non)symmetric
   matrix decompositions. Thus, SUPERLU must also be installed if you
   intend to use one of these classes. SUPERLU is available at this
   webpage (see also install-superlu.sh):

   http://crd-legacy.lbl.gov/~xiaoye/SuperLU/

4. UMFPACK:

   UMFPACK package can also be used to solve eigenvalue problems that
   require real or complex (non)symmetric/non-Hermitian matrix
   decompositions. UMFPACK is now part of the SuiteSparse package which
   can be obtained here (see also install-suitesparse.sh):

   http://faculty.cse.tamu.edu/davis/suitesparse.html

5. CHOLMOD

   CHOLMOD package is performing a Cholesky decomposition. Some of the
   symmetric problems can now interface with it. It is part of the
   SuiteSparse package which can be obtained here (see also
   install-suitesparse.sh):

   http://faculty.cse.tamu.edu/davis/suitesparse.html


## Documentation
## 
   Arpackpp user's manual is available in the doc directory. It contains
   all information needed to declare and solve eigenvalue problems using
   arpack++ classes and functions. Arpackpp computational modes and data
   types are also described in the manual. Instructions on how to
   install the above mentioned libraries are given in the INSTALL.md
   file. Moreover, README files were include in many arpackpp
   directories to give additional information about arpackpp files and
   examples.

## Using arpackpp:
## 
   As a collection of class templates, arpackpp need not to be compiled.
   Because templates are defined in header (.h) files, no object (.o) or
   library (.a) files have to be build, except those corresponding to
   other libraries required by arpackpp (see Dependencies above).
   Arpackpp header files are included in the "include" directory and can
   be moved to another directory if desired. An option in the form

   ```
   -I$(ARPACKPP_INC) \
   -I$(ARPACKPP_INC)/../examples/matrices \
   -I$(ARPACKPP_INC)/../examples/matprod
   ```

   should be added to the command line when compiling programs that use
   arpackpp. Here, ARPACKPP_INC is the name of the directory that
   contains all arpackpp header files. Note, depending on what type of
   problem you want so solve, you need to also include the example
   matrices and/or matprod directories (see examples). You can also use
   `cmake` (see below) with `make install` to install all headers to your
   system into a single directory.


## Compiling and running arpackpp examples:
## 
   Arpackpp supports cmake for the compilation of the examples. To build
   all examples, including the ones that depend on SuperLU, do

   ```
   $ mkdir ../arpackpp-build
   $ cd ../arpackpp-build
   $ cmake ../arpackpp -D SUPERLU=ON
   $ make examples
   ```

   For this to work all dependencies need to be installed (either on the
   system or in the external subdirectory). See INSTALL.md for details.
   Regular Makefiles (in-source build) are also still supported.


## Acknowledgements
## 
ARPACK++ authors:

-  Francisco M. Gomes (chico AT ime.unicamp.br)

-  Danny Sorensen (LASTNAME AT caam.rice.edu)

arpackpp (2.0.0 and above) authors:

-  Martin Reuter (LASTNAME AT mit.edu)



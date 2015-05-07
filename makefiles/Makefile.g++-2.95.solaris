# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# This file contains some definitions used to compile arpack++ examples
# with the g++ compiler under Sun Solaris.


# Defining the machine.

PLAT         = solaris

# Defining the compiler.

CPP          = g++

# Defining ARPACK++ directories.
# ARPACKPP_INC is the directory that contains all arpack++ header files.
# SUPERLU_DIR and UMFPACK_DIR must be set to ARPACKPP_INC.

ARPACKPP_DIR = $(HOME)/arpack++
ARPACKPP_INC = $(ARPACKPP_DIR)/include
SUPERLU_DIR  = $(ARPACKPP_INC)
UMFPACK_DIR  = $(ARPACKPP_INC)

# Defining ARPACK, LAPACK, UMFPACK, SUPERLU, BLAS and FORTRAN libraries.
# See the arpack++ manual or the README file for directions on how to 
# obtain arpack, umfpack and SuperLU packages. 
# UMFPACK_LIB and SUPERLU_LIB must be declared only if umfpack and superlu 
# are going to be used. Some BLAS and LAPACK fortran routines are 
# distributed along with arpack fortran code, but the user should verify 
# if optimized versions of these libraries are available before installing 
# arpack. The fortran libraries described below are those required to link
# fortran and c++ code using gnu g++ and f77 compiler under linux.
# Other libraries should be defined if the user intends to compile
# arpack++ on another environment.

ARPACK_LIB   = $(HOME)/lib/arpack_$(PLAT).a
LAPACK_LIB   = /opt/LAPACK/lapack_$(PLAT).a
UMFPACK_LIB  = $(HOME)/lib/umfpack_$(PLAT).a
SUPERLU_LIB  = $(HOME)/lib/superlu_$(PLAT).a
BLAS_LIB     = /opt/LAPACK/blas_$(PLAT).a 
FORTRAN_LIBS = /opt/SUNWspro/SC4.0/lib/libF77.a \
               /opt/SUNWspro/SC4.0/lib/libM77.a \
               /opt/SUNWspro/SC4.0/lib/libsunmath.a

# Defining g++ flags and directories.

# CPP_WARNINGS = -fpermissive 
CPP_WARNINGS = -Wall -ansi -pedantic-errors
CPP_DEBUG    = -g
CPP_OPTIM    = -O
CPP_LIBS     = 
CPP_INC      = /usr/local/include/g++-3

CPP_FLAGS    = $(CPP_DEBUG) -D$(PLAT) -I$(ARPACKPP_INC) -I$(CPP_INC) \
               $(CPP_WARNINGS)

# Putting all libraries together.

ALL_LIBS     = $(CPP_LIBS) $(ARPACK_LIB) $(LAPACK_LIB) \
               $(BLAS_LIB) $(FORTRAN_LIBS) 

# defining paths.

vpath %.h  $(ARPACK_INC)


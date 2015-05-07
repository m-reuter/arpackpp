# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/umfpack/nonsym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/nonsym

# compiling and linking all examples.

all: unsymreg unsymshf unsymgre unsymgsh unsymgsc usvd

# compiling and linking each nonsymmetric problem.

unsymreg:       unsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o unsymreg unsymreg.o $(UMFPACK_LIB) $(ALL_LIBS)

unsymshf:	unsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o unsymshf unsymshf.o $(UMFPACK_LIB) $(ALL_LIBS)

unsymgre:	unsymgre.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o unsymgre unsymgre.o $(UMFPACK_LIB) $(ALL_LIBS)

unsymgsh:	unsymgsh.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o unsymgsh unsymgsh.o $(UMFPACK_LIB) $(ALL_LIBS)

unsymgsc:	unsymgsc.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o unsymgsc unsymgsc.o $(UMFPACK_LIB) $(ALL_LIBS)

# compiling and linking svd problems.

usvd:	usvd.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o usvd usvd.o $(UMFPACK_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core unsymreg unsymshf unsymgre unsymgsh unsymgsc usvd

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(UMFPACK_DIR) -c $<



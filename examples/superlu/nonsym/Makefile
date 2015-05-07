# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/superlu/nonsym directory makefile

# including other makefiles.

include ../../../Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/nonsym

# compiling and linking all examples.

all: lnsymreg lnsymshf lnsymgre lnsymgsh lnsymgsc lsvd lsvd2

# compiling and linking each nonsymmetric problem.

lnsymreg:       lnsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lnsymreg lnsymreg.o $(SUPERLU_LIB) $(ALL_LIBS)

lnsymshf:	lnsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lnsymshf lnsymshf.o $(SUPERLU_LIB) $(ALL_LIBS)

lnsymgre:	lnsymgre.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lnsymgre lnsymgre.o $(SUPERLU_LIB) $(ALL_LIBS)

lnsymgsh:	lnsymgsh.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lnsymgsh lnsymgsh.o $(SUPERLU_LIB) $(ALL_LIBS)

lnsymgsc:	lnsymgsc.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lnsymgsc lnsymgsc.o $(SUPERLU_LIB) $(ALL_LIBS)

# compiling and linking svd problems.

lsvd:	lsvd.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsvd lsvd.o $(SUPERLU_LIB) $(ALL_LIBS)

lsvd2:	lsvd2.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsvd2 lsvd2.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core lnsymreg lnsymshf lnsymgre lnsymgsh lnsymgsc lsvd lsvd2

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<


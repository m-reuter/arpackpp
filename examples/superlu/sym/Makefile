# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/superlu/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining CSCMAT directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: lsymreg lsymshf lsymgreg lsymgshf lsymgbkl lsymgcay

# compiling and linking each symmetric problem.

lsymreg:	lsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsymreg lsymreg.o $(SUPERLU_LIB) $(ALL_LIBS)

lsymshf:	lsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsymshf lsymshf.o $(SUPERLU_LIB) $(ALL_LIBS)

lsymgreg:	lsymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsymgreg lsymgreg.o $(SUPERLU_LIB) $(ALL_LIBS)

lsymgshf:	lsymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsymgshf lsymgshf.o $(SUPERLU_LIB) $(ALL_LIBS)

lsymgbkl:	lsymgbkl.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsymgbkl lsymgbkl.o $(SUPERLU_LIB) $(ALL_LIBS)

lsymgcay:	lsymgcay.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o lsymgcay lsymgcay.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core lsymgbkl lsymgcay lsymgreg lsymgshf lsymreg lsymshf

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<



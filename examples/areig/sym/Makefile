# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/areig/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining areig directory.

AREIG_DIR = $(ARPACKPP_DIR)/examples/areig

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: asymreg asymshf asymgreg asymgshf asymgbkl asymgcay

# compiling and linking each symmetric problem.

asymreg:	asymreg.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o asymreg asymreg.o $(SUPERLU_LIB) $(ALL_LIBS)

asymshf:	asymshf.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o asymshf asymshf.o $(SUPERLU_LIB) $(ALL_LIBS)

asymgreg:	asymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o asymgreg asymgreg.o $(SUPERLU_LIB) $(ALL_LIBS)

asymgshf:	asymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o asymgshf asymgshf.o $(SUPERLU_LIB) $(ALL_LIBS)

asymgbkl:	asymgbkl.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o asymgbkl asymgbkl.o $(SUPERLU_LIB) $(ALL_LIBS)

asymgcay:	asymgcay.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o asymgcay asymgcay.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core asymgbkl asymgcay asymgreg asymgshf asymreg asymshf

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<


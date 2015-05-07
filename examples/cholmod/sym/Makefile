# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/cholmod/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining CSCMAT directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: csymreg csymshf csymgreg csymgshf #csymgbkl csymgcay

# compiling and linking each symmetric problem.

csymreg:	csymreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o csymreg csymreg.o $(CHOLMOD_LIB) $(ALL_LIBS)

csymshf:	csymshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o csymshf csymshf.o $(CHOLMOD_LIB) $(ALL_LIBS)

csymgreg:	csymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o csymgreg csymgreg.o $(CHOLMOD_LIB) $(ALL_LIBS)

csymgshf:	csymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o csymgshf csymgshf.o $(CHOLMOD_LIB) $(ALL_LIBS)

csymgbkl:	csymgbkl.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o csymgbkl csymgbkl.o $(CHOLMOD_LIB) $(ALL_LIBS)

csymgcay:	csymgcay.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o csymgcay csymgcay.o $(CHOLMOD_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core csymreg csymshf csymgreg csymgshf csymgbkl csymgcay

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) $(CHOLMOD_INC) -c $<



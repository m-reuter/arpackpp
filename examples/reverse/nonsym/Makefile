# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/reverse/nonsym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining objects.

EXMP_INC = $(ARPACKPP_DIR)/examples/matprod 
EXN_INC  = $(ARPACKPP_DIR)/examples/matprod/nonsym 

# compiling and linking all examples.

all: rnsymreg rnsymshf rnsymgre rnsymgsh rnsymgsc rsvd

# compiling and linking each nonsymmetric problem.

rnsymreg:       rnsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o rnsymreg rnsymreg.o \
            $(ALL_LIBS)

rnsymshf:       rnsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o rnsymshf rnsymshf.o \
            $(ALL_LIBS)

rnsymgre:       rnsymgre.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o rnsymgre rnsymgre.o \
            $(ALL_LIBS)

rnsymgsh:       rnsymgsh.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o rnsymgsh rnsymgsh.o \
            $(ALL_LIBS)

rnsymgsc:       rnsymgsc.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o rnsymgsc rnsymgsc.o \
            $(ALL_LIBS)

# compiling and linking svd problems.

rsvd:	rsvd.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o rsvd rsvd.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *.o *~ core rnsymreg rnsymshf rnsymgre rnsymgsh rnsymgsc rsvd

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -c $<


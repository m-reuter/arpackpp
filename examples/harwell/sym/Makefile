# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/harwell/sym directory makefile

# including other makefiles.

include ../../../Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: hsymstd hsymgen

# compiling and linking each symmetric problem.

hsymstd:	hsymstd.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o hsymstd hsymstd.o $(SUPERLU_LIB) $(ALL_LIBS)

hsymgen:	hsymgen.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o hsymgen hsymgen.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core hsymstd hsymgen

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<


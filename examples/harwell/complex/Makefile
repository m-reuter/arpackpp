# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/harwell/complex directory makefile

# including other makefiles.

include ../../../Makefile.inc

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/complex

# compiling and linking all examples.

all: hcompstd hcompgen

# compiling and linking each complex problem.

hcompstd:	hcompstd.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o hcompstd hcompstd.o $(SUPERLU_LIB) $(ALL_LIBS)

hcompgen:	hcompgen.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o hcompgen hcompgen.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core hcompstd hcompgen

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<


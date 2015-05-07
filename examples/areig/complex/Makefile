# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/areig/complex directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining areig directory.

AREIG_DIR = $(ARPACKPP_DIR)/examples/areig

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/complex

# compiling and linking all examples.

all: acompreg acompshf acompgre acompgsh

# compiling and linking each complex problem.

acompreg:        acompreg.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o acompreg acompreg.o $(SUPERLU_LIB) $(ALL_LIBS)

acompshf:        acompshf.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o acompshf acompshf.o $(SUPERLU_LIB) $(ALL_LIBS)

acompgsh:        acompgsh.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o acompgsh acompgsh.o $(SUPERLU_LIB) $(ALL_LIBS)

acompgre:        acompgre.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o acompgre acompgre.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core acompgre acompgsh acompreg acompshf

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -c $<


# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/dense/complex directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining DENSEMAT directory.

DENSEMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/complex

# compiling and linking all examples.

all: dcompreg dcompshf dcompgre dcompgsh

# compiling and linking each complex problem.

dcompreg:       dcompreg.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dcompreg dcompreg.o $(ALL_LIBS)

dcompshf:	dcompshf.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dcompshf dcompshf.o $(ALL_LIBS)

dcompgre:	dcompgre.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dcompgre dcompgre.o $(ALL_LIBS)

dcompgsh:	dcompgsh.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dcompgsh dcompgsh.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core dcompreg dcompshf dcompgre dcompgsh

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -c $<



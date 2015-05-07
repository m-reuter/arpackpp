# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/dense/nonsym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining densemat directory.

DENSEMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/nonsym

# compiling and linking all examples.

all: dnsymreg dnsymshf dnsymgre dnsymgsh dnsymgsc dsvd

# compiling and linking each nonsymmetric problem.

dnsymreg:       dnsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dnsymreg dnsymreg.o $(ALL_LIBS)

dnsymshf:	dnsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dnsymshf dnsymshf.o $(ALL_LIBS)

dnsymgre:	dnsymgre.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dnsymgre dnsymgre.o $(ALL_LIBS)

dnsymgsh:	dnsymgsh.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dnsymgsh dnsymgsh.o $(ALL_LIBS)

dnsymgsc:	dnsymgsc.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dnsymgsc dnsymgsc.o $(ALL_LIBS)

dsvd:	dsvd.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsvd dsvd.o $(ALL_LIBS)

dsvd2:	dsvd2.o
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -o dsvd2 dsvd2.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core dnsymreg dnsymshf dnsymgre dnsymgsh dnsymgsc dsvd

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(DENSEMAT_DIR) -c $<

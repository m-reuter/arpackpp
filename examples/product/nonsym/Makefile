# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/product/nonsym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining objects.

EXMP_INC = $(ARPACKPP_DIR)/examples/matprod 
EXN_INC  = $(ARPACKPP_DIR)/examples/matprod/nonsym 

# compiling and linking all examples.

all: nsymreg nsymshf nsymgreg nsymgshf nsymgscr nsymgsci svd

# compiling and linking each nonsymmetric problem.

nsymreg: 	nsymreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o nsymreg nsymreg.o $(ALL_LIBS)

nsymshf:	nsymshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o nsymshf nsymshf.o $(ALL_LIBS)

nsymgreg:	nsymgreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o nsymgreg nsymgreg.o $(ALL_LIBS)

nsymgshf:	nsymgshf.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o nsymgshf nsymgshf.o $(ALL_LIBS)

nsymgscr:	$(COMM_H) nsymgscr.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o nsymgscr nsymgscr.o $(ALL_LIBS)

nsymgsci:	$(COMM_H) nsymgsci.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o nsymgsci nsymgsci.o $(ALL_LIBS)

# compiling and linking svd problems.

svd:	svd.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -o svd svd.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core nsymgreg nsymgsci nsymgscr nsymgshf nsymreg nsymshf svd

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXN_INC) -c $<


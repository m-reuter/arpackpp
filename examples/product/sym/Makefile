# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/product/sym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

# defining objects.
 
EXMP_INC = $(ARPACKPP_DIR)/examples/matprod
EXS_INC  = $(ARPACKPP_DIR)/examples/matprod/sym

# compiling and linking all examples.

all: symreg symshft symgreg symgshft symgbklg symgcayl

# compiling and linking each symmetric problem.

symreg:		symreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o symreg symreg.o $(ALL_LIBS) 

symshft:	symshft.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o symshft symshft.o $(ALL_LIBS)

symgreg:	symgreg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o symgreg symgreg.o $(ALL_LIBS)

symgshft:	symgshft.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o symgshft symgshft.o $(ALL_LIBS)

symgbklg:	symgbklg.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o symgbklg symgbklg.o $(ALL_LIBS)

symgcayl:	symgcayl.o
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -o symgcayl symgcayl.o $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core symgbklg symgcayl symgreg symgshft symreg symshft

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(EXMP_INC) -I$(EXS_INC) -c $<


FOR = f95 #/usr/local/intel/Compiler/11.1/064/bin/intel64/ifort 
CC  = g++
#note use -mp -prec_div instead of -fast for accurate math calculations
OPTFLAGS = -O3 
PRECFLAG =  #-r16
#OPTFLAGS =  -g -C -debug variable_locations   -traceback  -fpe0
RM  = rm -f
FPP = -cpp
INCLUDE  = -I/usr/include/eigen3/
FPPFLAGS = #-DVERBOSE #-DWRITE_QINFO # -DDEBUG 
#OPEN_MP   = -openmp

FORFLAGS = $(PRECFLAG) $(OPTFLAGS) $(FPP) $(FPPFLAGS) $(OPEN_MP) $(INCLUDE)  
CCFLAGS = $(OPTFLAGS) $(INCLUDE)

LIBS  = -llapack -lblas -lgfortran -ljsoncpp


.SUFFIXES:
 
.SUFFIXES:	.o .f90 .f .c .cpp
.f90.o:
	$(FOR) $(FORFLAGS) -c  $<
.f.o:
	$(FOR) $(FORFLAGS) -c  $<
.c.o:
	$(CC) $(CCFLAGS) -c $<
.cpp.o:
	$(CC) $(CCFLAGS) -c $<


INC = Params.hpp

TAGT = ising

OBJ0 = 
OBJ1 = 
OBJ2 =
OBJ3 = 
OBJ4 = 
OBJ5 = 
OBJ6 = 
OBJ7 = 
OBJ8 = 
OBJ9 = ising.o

PACK = ising

all: $(TAGT) 

$(TAGT) : $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ4) $(OBJ5) $(OBJ6) \
          $(OBJ7) $(OBJ8) $(OBJ9) $(INC) makefile
	$(CC) -o $(TAGT) $(FORFLAGS) $(OBJ0) $(OBJ1) $(OBJ2) \
        $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(LIBS) 

$(OBJ0): $(INC)
$(OBJ1): $(OBJ0) $(INC)
$(OBJ2): $(OBJ1) $(OBJ0) $(INC)
$(OBJ3): $(OBJ2) $(OBJ1) $(OBJ0) $(INC)
$(OBJ4): $(OBJ3) $(OBJ2) $(OBJ1) $(OBJ0) $(INC)
$(OBJ5): $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1) $(OBJ0) $(INC)
$(OBJ6): $(OBJ5) $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1) $(OBJ0) $(INC)
$(OBJ7): $(OBJ6) $(OBJ5) $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1) $(OBJ0) $(INC)
$(OBJ8): $(OBJ7) $(OBJ6) $(OBJ5) $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1) $(OBJ0) $(INC)
$(OBJ9): $(OBJ8) $(OBJ7) $(OBJ6) $(OBJ5) $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1) $(OBJ0) $(INC)


clean :
	$(RM) *.mod *.o  $(TAGT)


tarfile:
	tar cvfhz $(PACK)_`hostname`_`date '+%d%b%y-%p'`.tgz \
        makefile $(INC) zn.cpp in.json


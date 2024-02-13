#!/bin/sh
 CC = gfortran
 #FFLAGS= -g -fbacktrace -ffpe-trap=zero,overflow,underflow -c -O3  
ODIR=src/objs
FFLAGS= -c -O3 -fno-range-check
FFLAGS2=-L/usr/lib/x86_64-linux-gnu/lapack/curren -llapack -L/usr/lib/x86_64-linux-gnu/blas/current -lblas
#FFLAGS = -c -O0 -g -Wall -Wextra -pedantic -fbounds-check -fno-range-check -fbacktrace
#FFLAGS2= -O3 -g -Wall -Wextra -pedantic -fbounds-check
#FFLAGS2= -O3 
#FFLAGS2= -pg #for gprof
_OBJS= types.o readMakeInp.o math.o dsyev.o jacobiMethod.o writting.o main.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
gauss_jordan :  ${OBJS} 
	$(CC) ${OBJS}  -o   ./huckelTri ${FFLAGS2}

$(ODIR)/%.o : src/%.f90
	$(CC) ${FFLAGS} $< -o  $@  -J$(ODIR)

clean: 
	rm $(ODIR)/*.o 
	rm $(ODIR)/*.mod




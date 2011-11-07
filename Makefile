.SUFFIXES : .o .c .f .f90
#
FC = gfortran
CC = cc

# optimization flags
OFLAG = -O3

# Fortran and C flags
FFLAGS = $(DFLAGS) $(OFLAG) $(DEBUG)
CFLAGS = $(DFLAGS) $(OFLAG) $(DEBUG)
 
# Fortran 90/95 format
FREE = -ffree-form -ffree-line-length-none
 
# Debug flags
#DEBUG  = -g -traceback #-D_DEBUG
#
# MPI or OpenMP
#MPI = -DOMP -fopenmp
#
BASE = strfac
MAIN = ${BASE}
#
# source and rules
#====================================================================
SOURCE = prec.o error.o  vardef.o identify.o selection.o ha.o       \
         strfac.o separate.o sfoneq.o paircorr.o pppc.o output.o    \
         csro.o sisf.o strfacdist.o main.o

all:  ${MAIN}

${MAIN}:  $(SOURCE)
	$(FC) $(FREE) $(SOURCE) $(OFLAG) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${MAIN}

tar:
	rm -f ${BASE}.tar; tar -czvf ${BASE}.tar.gz *.f90 Makefile MAKE/Make.* README


.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<

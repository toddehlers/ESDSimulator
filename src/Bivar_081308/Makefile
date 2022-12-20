.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

# using GNU compiler:
FCOMP = gfortran
CCOMP = gcc

# build with all cheks, does give a lot of errors!
#PGFLAGS = -c -O3 -Wall -Werror -Warray-bounds -Wconversion -Wunderflow -Wsurprising -fbacktrace -fdump-core

# normal build:
# PGFLAGS = -c -O3 -ffree-line-length-none -fopenmp
PGFLAGS = -c -O3 -ffree-line-length-none -Wall -cpp -fno-automatic -ffpe-trap=invalid,overflow
CFLAGS = -c -O3 -Wall

# debug build
#PGFLAGS = -c -g3 -ffree-line-length-none -fopenmp -ffpe-trap=invalid,zero,overflow,underflow
#CFLAGS = -c -g3


# profile build: option -pg for compiler and linker
# visualize profile:
# gprof pecube | gprof2dot.py | dot -Tpng -o profile.png
# http://code.google.com/p/jrfonseca/wiki/Gprof2Dot


# using intel compiler
#FCOMP = ifort
# optimize for newer processors
#PGFLAGS = -c -fast
# for all other i386 processors:
#PGFLAGS = -c -O2


OBJECTS = \
bivar.o \
bivariate_interpolation.o


# needed because the make command thinks .mod files are compiled with a MODULA compiler
%.o : %.mod

.f.o:
	$(FCOMP) $(PGFLAGS) $*.f

.f90.o:
	$(FCOMP) $(PGFLAGS) $*.f90

.c.o:
	$(CCOMP) $(CFLAGS) $*.c

all: bivar_linux

bivar_linux: $(OBJECTS)
	$(FCOMP) -O3  *.o -o ../../bivar_linux
#	$(FCOMP) -O3 -fopenmp *.o -o ../../bivar_linux

clean:
	rm -f *.o *.mod


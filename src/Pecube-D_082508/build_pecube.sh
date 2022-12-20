#!/bin/bash
#
# build_pecube.sh - compiles pecube executable
#
# Portland Group Fortran compiler
#pgf90 -O3 -c pecube_src_files/nrtype.f90
#pgf90 -O3 -c pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#pgf90 -O3 -c pecube_src_files/*.f90 pecube_src_files/*.f
#pgf90 -O3 -o pecube *.o

# Portland Group Fortran compiler safe
#pgf90 -c pecube_src_files/nrtype.f90
#pgf90 -c pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#pgf90 -c pecube_src_files/*.f90 pecube_src_files/*.f
#pgf90 -o pecube *.o

# Portland Group Fortran compiler w/ debug/warnings
#pgf90 -C -w -c pecube_src_files/nrtype.f90
#pgf90 -C -w -c pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#pgf90 -C -w -c pecube_src_files/*.f90 pecube_src_files/*.f
#pgf90 -C -w -o pecube *.o

# Intel Fortran compiler
#ifort -O3 -c nrtype.f90
#ifort -O3 -c nr.f90 nrutil.f90
#ifort -O3 -c *.f90 *.f
#ifort -O3 -o pecube *.o

# Intel Fortran compiler optimized for Core 2 Duos
#ifort -fast -c pecube_src_files/nrtype.f90
#ifort -fast -c pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#ifort -fast -c pecube_src_files/*.f90 pecube_src_files/*.f
#ifort -fast -o pecube *.o

# Intel Fortran compiler w/ debug/warning flags
#ifort -c -C -w -debug all pecube_src_files/nrtype.f90
#ifort -c -C -w -debug all pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#ifort -c -C -w -debug all pecube_src_files/*.f90 pecube_src_files/*.f
#ifort -o pecube -C -w -debug all *.o

# Intel Fortran compiler for idb debugger
#ifort -c -g -pg pecube_src_files/nrtype.f90
#ifort -c -g -pg pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#ifort -c -g -pg pecube_src_files/*.f90 pecube_src_files/*.f
#ifort -o pecube -g -pg *.o

# Latest GNU Fortran compiler (slower than ifort or pgf90, but free)
#gfortran -O3 -Wall -Werror -Warray-bounds -Wconversion -Wunderflow -Wsurprising -fbacktrace -fdump-core -c nrtype.f90
#gfortran -O3 -Wall -Werror -Warray-bounds -Wconversion -Wunderflow -Wsurprising -fbacktrace -fdump-core -c nr.f90 nrutil.f90
#gfortran -O3 -Wall -Werror -Warray-bounds -Wconversion -Wunderflow -Wsurprising -fbacktrace -fdump-core -c *.f90 *.f
#gfortran -O3 -Wall -Werror -Warray-bounds -Wconversion -Wunderflow -Wsurprising -fbacktrace -fdump-core -o pecube *.o

gfortran -O3 -ffree-line-length-none -c nrtype.f90
gfortran -O3 -ffree-line-length-none -c nr.f90 nrutil.f90
gfortran -O3 -ffree-line-length-none -c *.f90 *.f
gfortran -O3 -ffree-line-length-none -o pecube *.o

# Latest GNU Fortran compiler w/ debug/warnings
#gfortran -C -w -c pecube_src_files/nrtype.f90
#gfortran -C -w -c pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#gfortran -C -w -c pecube_src_files/*.f90 pecube_src_files/*.f
#gfortran -C -w -o pecube *.o

# Another open source Fortran compiler
#g95 -c pecube_src_files/nrtype.f90
#g95 -c pecube_src_files/nr.f90 pecube_src_files/nrutil.f90
#g95 -c pecube_src_files/*.f90 pecube_src_files/*.f
#g95 -o pecube *.o

file='../../input'
if [ -e $file ]
  then cp -i pecube ../../
fi

# Remove object files after compilation
rm -f *.o
rm -f *.mod

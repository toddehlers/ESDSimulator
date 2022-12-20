#Makefile for program cascade
# developed by Jean Braun
#              Research School of Earth Sciences
#              The Australian National University
#              Canberra, ACT 0200, Australia
#              Tel:+61-2-6249-5512
#              Fax:+61-2-6249-5443
#              email: Jean.Braun@anu.edu.au

# compilation flags (edit here to optimize for your machine...)
# the -p option turns the profiling option on

#FCOMP = gfortran
#FCOMP = g77
FCOMP = ifort
#FCOMP = pgf77
CCOMP = gcc
#CCOMP = icc

#PGFLAGS = -c -Mprof=func -O2
#CFLAGS = -c -O2
#LFLAGS = -Mprof=func

# Optimization level 3
PGFLAGS = -c -O3
CFLAGS = -c -O3
LFLAGS = 

# Fast Intel flags for C2D machines
#PGFLAGS = -c -fast
#CFLAGS = -c -fast
#LFLAGS = -fast

# Debugging flags
#PGFLAGS = -c -w -C
#CFLAGS = -c -w -C
#CFLAGS = -c
#LFLAGS =

# this where the X11 include files must be (system dependent)
# if your system is properly setup, leave it blank

INCLUDE = /usr/openwin/include

# object (all the interesting bits to do the erosion/sedimentation
#         computations)

OBJECTS = \
check_for_removal.o \
change_sea_level.o \
fluvial_erosion.o \
find_neighbours.o \
find_donors.o \
find_order.o \
random.o \
find_surface.o \
check_mesh.o \
debug.o \
diffusion_erosion.o \
solve_diffusion.o \
build_a.o \
erosional_properties.o \
tectonic_uplift.o \
find_catchment.o \
flexure.o \
mdian2.o \
read_but_skip_comment.o \
iread_but_skip_comment.o \
initialize_nodal_geometry.o \
initialize_general_parameters.o \
tectonic_movement.o \
update_bedrock.o \
update_time_step.o \
write_tecplot_output.o \
write_output.o \
rainmaker.o \
show.o	\
landslide.o \
landslide_simple.o \
find_dslope.o

# OBJECTS_XPLOT = \
# LIB/XPLOT/xplot.o \
# LIB/XPLOT/psplot.o
# alternate dummy graphic routines
OBJECTS_XPLOT = \
LIB/XPLOT_ALT/xplot.o \
LIB/XPLOT_ALT/psplot.o

ARCHIVE_XPLOT = \
xplot.o \
psplot.o

OBJECTS_NN2D = \
LIB/NN2D/del_flip.o \
LIB/NN2D/del_sub.o \
LIB/NN2D/delaun.o \
LIB/NN2D/nn1.o \
LIB/NN2D/nn2.o \
LIB/NN2D/nn_remove.o \
LIB/NN2D/nnplot.o \
LIB/NN2D/qhullf_dummy.o \
LIB/NN2D/stack.o \
LIB/NN2D/stackpair.o \
LIB/NN2D/volume.o

ARCHIVE_NN2D = \
del_flip.o \
del_sub.o \
delaun.o \
nn1.o \
nn2.o \
nn_remove.o \
nnplot.o \
qhullf_dummy.o \
stack.o \
stackpair.o \
volume.o

# all the non-interesting bits

UTILS = \
sinft.o \
four1.o \
realft.o

# compilation only

.f.o:
	$(FCOMP) $(PGFLAGS) $*.f
.c.o:
	$(CCOMP) $(CFLAGS) $*.c

# create an executable
# note that the option -lsocket is something that we have found necessary
# on our Sun workstations; it is usually not needed on a well maintained
# Unix system...
cascade: $(OBJECTS) $(UTILS) cascade.f cascade.h
	$(FCOMP) cascade.f $(LFLAGS) $(OBJECTS) $(UTILS) -LLIB -lnn2d -lxplot -o cascade
	if [ -e '../../input' ]; then cp -i cascade ../../; fi


lib: $(OBJECTS_XPLOT) $(OBJECTS_NN2D)
	ar rcv LIB/libnn2d.a $(ARCHIVE_NN2D)
	ar rcv LIB/libxplot.a $(ARCHIVE_XPLOT)

clean:
	rm -f *.o
	rm -f LIB/*.a

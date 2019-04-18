.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .e .r .y .yr .ye .l .p .sh .csh .h

#--- COMPILER CONFIGURATION ----------------------------------------------------
# Only one set of these variables should be defined in this file. Only one set
# of lines should start with F90FLAGS, F77FLAGS, etc. down to the line that
# reads "OBJECTS = \". All others should start with a #. Remove the # to enable
# a given set of compilers.

#--- PECUBE FOR FORWARD MODELLING ----------------------------------------------
# The default behavior is for Pecube to build for running forward models. If you
# are interested in using Pecube for inverse modelling, please comment out the
# compiler options here and try one of the options in the section below.

# Default compiler options (slow, choose a better set below)
F90FLAGS = -c -O0
F77FLAGS = -c -O0 -fd-lines-as-comments -ffixed-line-length-132
CFLAGS = -c -O0
LDFLAGS = 
FF90 = gfortran
FF77 = gfortran
CC = gcc

# Codecov compiler options (slow, choose a better set below)
#F90FLAGS = -c -O0 -fprofile-arcs -ftest-coverage
#F77FLAGS = -c -O0 -fd-lines-as-comments -ffixed-line-length-132 -fprofile-arcs -ftest-coverage
#CFLAGS = -c -O0 -fprofile-arcs -ftest-coverage
#LDFLAGS = -lgcov --coverage
#FF90 = gfortran
#FF77 = gfortran
#CC = gcc

# GNU compilers (e.g., gcc, gfortran)
#F90FLAGS = -c -O2
#F77FLAGS = -c -O2 -fd-lines-as-comments -ffixed-line-length-132
#CFLAGS = -c -O2
#LDFLAGS =
#FF90 = /usr/local/bin/gfortran
#FF77 = /usr/local/bin/gfortran
#CC = /usr/local/bin/gcc

# Intel compilers (e.g., icc, ifort))
#F90FLAGS = -c -O2
#F77FLAGS = -c -O2
#CFLAGS = -c -O2
#LDFLAGS =
#FF90 = /usr/local/intel/composerxe-2011.0.084/bin/intel64/ifort
#FF77 = /usr/local/intel/composerxe-2011.0.084/bin/intel64/ifort
#CC = /usr/local/intel/composerxe-2011.0.084/bin/intel64/icc

#--- PECUBE FOR INVERSE MODELLING ----------------------------------------------
# NOTE: In order to compile Pecube for inverse modelling you need to replace the
# na.f file with another version. The options are na.f.ORIG for the original NA
# code or na.f.newNA for the more efficient updated version (recommended). To
# use one of these options, simply rename the original na.f file and replace it
# with one of the options. For example:
#
# mv na.f na.f.serial
# mv na.f.newNA na.f

# Parallel version with GNU compilers
#F90FLAGS = -c -O2
#F77FLAGS = -c -O2 -cpp -DNA_MPIP=1 -DNA_F77_COMPATIBILITY=0
#CFLAGS = -c -O2
#FF90 = /usr/local/openmpi.gcc-1.2.9/bin/mpif90
#FF77 = /usr/local/openmpi.gcc-1.2.9/bin/mpif77
#CC = /usr/local/openmpi.gcc-1.2.9/bin/mpicc

# Parallel version with Intel compilers
#F90FLAGS = -c -O2
#F77FLAGS = -c -O2 -fpp -DNA_MPIP=1 -DNA_F77_COMPATIBILITY=0
#CFLAGS = -c -O2
#FF90 = /usr/local/openmpi.intel/bin/mpif90
#FF77 = /usr/local/openmpi.intel/bin/mpif77
#CC = /usr/local/openmpi.intel/bin/mpicc

OBJECTS = \
Mad_He.o \
Mad_Trax.o \
Mad_Trax_Zircon.o \
Pecube.o \
calc_raman.o \
calculate_fault_parameters.o \
calculate_topo_advection.o \
calculate_ages.o \
create_pecube_in.o \
define_proc.o \
find_distance.o \
find_dt.o \
find_element.o \
find_neighbours.o \
find_temp.o \
find_triangle.o \
find_velo.o \
find_velocity.o \
four1.o \
interpol3d.o \
interpolate.o \
isostatic_rebound.o \
make_matrix.o \
move_fault.o \
na.o \
na_utils.o \
read_in_fault_parameters.o \
realft.o \
screen_counter.o \
sinft.o \
solve_iterative.o \
tridag.o \
KETCH.o

OBJECTS_TEST = \
testFault.o \
calculate_fault_parameters.o \
calculate_topo_advection.o \
find_distance.o \
find_velocity.o \
move_fault.o \
read_in_fault_parameters.o \
testTopo.o

OBJECTS_VTK = \
vtk.o

OBJECTS_NA = \
showNA.o \
create_pecube_in.o \
read_in_fault_parameters.o \
calculate_fault_parameters.o \
calculate_topo_advection.o \
find_triangle.o

OBJECTS_CASCADE = \
cascade2Pecube.o

.f90.o:
	$(FF90) $(F90FLAGS) $*.f90

.F90.o:
	$(FF90) $(F90FLAGS) $*.F90

.f.o:
	$(FF77) $(F77FLAGS) $*.f

.c.o:
	$(CC) $(CFLAGS) $*.c

Pecube:	$(OBJECTS)
	$(FF90) -O2 $(OBJECTS) -o ../bin/Pecube $(LDFLAGS)

Test: $(OBJECTS_TEST)
	$(FF90) -O2 $(OBJECTS_TEST) -o ../bin/Test $(LDFLAGS)

Vtk: $(OBJECTS_VTK)
	$(FF90) -O2 $(OBJECTS_VTK) -o ../bin/Vtk $(LDFLAGS)

ShowNA:	$(OBJECTS_NA)
	$(FF90) -O2 $(OBJECTS_NA) -o ../bin/ShowNA $(LDFLAGS)

Cascade2Pecube:	$(OBJECTS_CASCADE)
	$(FF90) -O2 $(OBJECTS_CASCADE) -o ../bin/Cascade2Pecube $(LDFLAGS)

all: $(OBJECTS) $(OBJECTS_TEST)  $(OBJECTS_VTK)
	make Pecube
	make Test
	make Vtk
	make ShowNA
	make Cascade2Pecube

clean:
	rm *.o *.mod

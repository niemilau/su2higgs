# parallel makefile, requires MPI

CC := mpicc -O3

#CC := mpicc -ggdb3 -g -O3 -D MPI

## CFLAGS is implicitly added after OBJECTS. Use -D flags to control what fields / features are included:
# -DU1 : include compact hypercharge field
# -DNHIGGS=<int> : number of Higgs doublets, but only 0, 1, 2 are supported
# -DTRIPLET : include adjoing Higgs field
# -DSINGLET : include real gauge singlet scalar
# -DCORRELATORS : measure some two-point functions
# -DBLOCKING : do blocking transformations on the lattice to reduce noise (with correlation measurements only)
# -DGRADFLOW : do gradient flow smoothing
#
# Note that not all of the above flags work together.
 
#CFLAGS := -D MPI -D HIGGS -D TRIPLET #-D U1
#CFLAGS := -D MPI -D HIGGS -D TRIPLET
CFLAGS := -D MPI -D NHIGGS=1 -D SINGLET -D U1
#CFLAGS := -D MPI -D TRIPLET -D CORRELATORS -D BLOCKING -D GRADFLOW

LIBS := -lm

OBJECTS := src/main.o src/generic/mersenne.o src/layout.o src/comms.o src/alloc.o src/init.o src/parameters.o src/su2u1.o src/staples.o src/measure.o \
	src/update.o src/checkpoint.o src/metropolis.o src/heatbath.o src/overrelax.o src/multicanonical.o \
	src/blocking.o src/z_coord.o src/magfield.o src/gradflow.o src/correlation.o src/hb_trajectory.o

BINARY := build/su2

FILES := $(OBJECTS) $(BINARY)

all: su2

su2: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm $(OBJECTS)

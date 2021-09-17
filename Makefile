# parallel makefile, requires MPI

CC := mpicc -g -O3 -march=native

#CC := mpicc -ggdb3 -g -O3 -D MPI

## CFLAGS is implicitly added after OBJECTS, so make sure to use the correct CFLAGS...
 
#CFLAGS := -D MPI -D HIGGS -D TRIPLET #-D U1
#CFLAGS := -D MPI -D HIGGS -D TRIPLET
CFLAGS := -D MPI -D NHIGGS=1 -D U1 -D SINGLET -D CORRELATORS -D BLOCKING
#CFLAGS := -D MPI -D TRIPLET -D CORRELATORS -D BLOCKING -D GRADFLOW

LIBS := -lm

OBJECTS := main.o generic/mersenne.o layout.o comms.o alloc.o init.o parameters.o su2u1.o staples.o measure.o \
	update.o checkpoint.o metropolis.o heatbath.o overrelax.o multicanonical.o \
	smearing.o blocking.o z_coord.o magfield.o gradflow.o correlation.o hb_trajectory.o

BINARY := build/su2

FILES := $(OBJECTS) $(BINARY)

all: su2

su2: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm $(OBJECTS)

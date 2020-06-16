# parallel makefile, requires MPI

CC := mpicc -g -O3 -march=native

#CC := mpicc -ggdb3 -g -O3 -D MPI

## CFLAGS is implicitly added after OBJECTS, so make sure to use the correct CFLAGS...
 
#CFLAGS := -D MPI -D HIGGS -D TRIPLET #-D U1
#CFLAGS := -D MPI -D HIGGS -D TRIPLET
CFLAGS := -D MPI -D HIGGS -D TRIPLET -D GRADFLOW -D MEASURE_Z

LIBS := -lm

OBJECTS := main.o layout.o comms.o alloc.o init.o parameters.o su2u1.o measure.o \
	update.o checkpoint.o metropolis.o heatbath.o overrelax.o multicanonical.o \
	z_coord.o magfield.o gradflow.o

BINARY := build/su2

FILES := $(OBJECTS) $(BINARY)

all: su2

su2: $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm *.o

# parallel makefile, requires MPI

#CC := gcc -g -O3 -D HIGGS #-D TRIPLET #-D DEBUG

CC := mpicc -g -O3 -march=native

#CC := mpicc -ggdb3 -g -O3 -D MPI

## CFLAGS is implicitly added after OBJECTS, so make sure to use the correct CFLAGS...
 
#CFLAGS := -D MPI -D HIGGS -D TRIPLET #-D U1
#CFLAGS := -D MPI -D HIGGS -D TRIPLET -D WALL
CFLAGS := -D MPI -D HIGGS -D TRIPLET -D WALL

LIBS := -lm

OBJECTS := main.o layout.o comms.o alloc.o init.o parameters.o su2u1.o measure.o \
	update.o checkpoint.o metropolis.o heatbath.o overrelax.o multicanonical.o magfield.o

OBJECTS_WALL := wallprofile.o

BINARY := build/su2

FILES := $(OBJECTS) $(BINARY)

all: su2

su2: $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

su2wall: $(OBJECTS) $(OBJECTS_WALL)
	$(CC) $(CFLAGS) $(OBJECTS) $(OBJECTS_WALL) $(LIBS) -o build/su2wall

clean:
	rm *.o

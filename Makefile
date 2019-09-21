# parallel makefile, requires MPI

#CC := gcc -g -O3 -D HIGGS #-D TRIPLET #-D DEBUG

CC := mpicc -g -O3 -D MPI -D HIGGS -D TRIPLET -D U1

#CC := mpicc -ggdb3 -g -O3 -D MPI

CFLAGS := 

LIBS := -lm

OBJECTS := main.o layout.o comms.o alloc.o init.o parameters.o su2u1.o measure.o \
	update.o checkpoint.o metropolis.o heatbath.o overrelax.o multicanonical.o

BINARY := build/su2

FILES := $(OBJECTS) $(BINARY)

all: su2

su2: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#

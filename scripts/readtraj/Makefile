# parallel makefile, requires MPI

CC := gcc -g -O3

CFLAGS :=

LIBS := -lm

OBJECTS := main.o readtraj.o

BINARY := readtraj

FILES := $(OBJECTS) $(BINARY)

all: readtraj

readtraj: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#

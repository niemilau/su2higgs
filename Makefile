# parallel makefile, requires MPI

## Use -D flags to control what fields / features are included:
# -DU1 : include compact hypercharge field
# -DNHIGGS=<int> : number of Higgs doublets, but only 0, 1, 2 are supported
# -DTRIPLET : include adjoing Higgs field
# -DSINGLET : include real gauge singlet scalar
# -DCORRELATORS : measure some two-point functions
# -DBLOCKING : do blocking transformations on the lattice to reduce noise (with correlation measurements only)
# -DGRADFLOW : do gradient flow smoothing
#
# Note that not all of the above flags work together.

PROGRAM_CFLAGS := -DNHIGGS=1 -DSINGLET -DU1
#PROGRAM_CFLAGS := -DHIGGS -DTRIPLET
#PROGRAM_CFLAGS := -DTRIPLET -DCORRELATORS -DBLOCKING -DGRADFLOW


## For non-MPI builds use 'make SERIAL=1'
ifdef SERIAL
	CFLAGS := $(PROGRAM_CFLAGS)
	CC := gcc -O3
else
	CFLAGS := -DMPI $(PROGRAM_CFLAGS)
	CC := mpicc -O3
endif

LIBS := -lm


SRC_DIR := src
BUILD_DIR := build
BINARY_DIR := bin

SOURCES := main.c generic/mersenne.c layout.c comms.c alloc.c init.c parameters.c su2u1.c staples.c measure.c \
	update.c checkpoint.c metropolis.c heatbath.c overrelax.c multicanonical.c \
	blocking.c z_coord.c magfield.c gradflow.c correlation.c hb_trajectory.c

OBJECTS := $(addprefix $(BUILD_DIR)/,$(SOURCES:.c=.o))

BINARY := bin/su2


.PHONY: makedirs all clean

all: makedirs $(BINARY)

makedirs:
	mkdir -p $(BUILD_DIR) $(BUILD_DIR)/generic $(BINARY_DIR)

$(BINARY): $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

# Compile sources to .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

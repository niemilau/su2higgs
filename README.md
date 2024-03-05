## su2higgs

# Compiling

Simply run ```make``` in the repository root. The resulting binary will go to ```bin``` directory. By default we build a parallel MPI program (requires ```Open MPI```). A non-parallel, single-core program can be compiled with ```make SERIAL=1```.

```-D``` flags can be used to specify what field content and other features to include in simulations. See the Makefile for examples.

# Running

The program takes a configuration file as command line argument. This is where you can specify things like lattice size, input parameters to the lattice action (in ```a=1``` units) and whether a multicanonical algorithm should be used.  A sample ```config``` file is included in the repo.

Example run with 4 cores:
```mpirun -np 4 bin/su2 config```

The program computes volume averages (averages over all lattice sites) of local operators and stores them in plain text file ```measure``` (name can be changed in the config file). It also produces a ```labels``` file containing column labels for the measurement file. Individual field configurations are only stored at infrequent checkpoints that store a snapshot of the lattice system in a binary file (default name: ```lattice```).

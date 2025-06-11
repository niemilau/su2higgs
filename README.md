# su2higgs

This is an MPI-powered parallel C-program for Monte Carlo simulations of statistical and quantum field theories with SU(2) gauge symmetry and Higgs fields in different representations. Particular focus is placed on studying first-order phase transitions with multicanonical sampling.


## Supported field content

- SU(2) gauge fields with standard Wilson action
- Compact U(1) gauge fields
- SU(2) x U(1) fundamental Higgs
- SU(2) adjoint Higgs
- Real gauge singlet scalar

## Compiling

Simply run ```make``` in the repository root. The resulting binary will go to ```bin``` directory. By default we build a parallel MPI program (requires ```OpenMPI```). A non-MPI, serial program can be compiled with ```make SERIAL=1```.

```-D``` flags can be used to specify what field content and other features to include in simulations. See the Makefile for examples.

Compiles without warnings on GCC 9.4.0.

## Running

The program takes a configuration file as command line argument. This is where you can specify things like lattice size, input parameters to the lattice action (in ```a=1``` units) and whether a multicanonical algorithm should be used.  A sample config file is included in the repo.

Example run with 4 cores:
```mpirun -np 4 bin/su2 config```

The program computes volume averages (averages over all lattice sites) of local operators and stores them in plain text file 'measure' (name changeable in the config file). It also produces a 'labels' file containing column labels for the measurement file. Individual field configurations are only stored at infrequent checkpoints that store a snapshot of the lattice system in a binary file (default name: 'lattice').

## Literature

For physics background see https://arxiv.org/abs/2405.01191 and references therein.

This code has been used in the following scientific publications:
- Niemi et al, Phys.Rev.Lett. 126 (2021) 17, 171802 (arXiv: 2005.11332)
- Niemi et al, JHEP 02 (2023) 212 (arXiv: 2206.14487)
- Niemi et al, Phys.Rev.D 110 (2024) 11, 115016 (arXiv: 2405.01191)

If you use this code in your scientific work, please consider citing any of the above and including a link to this Github page.

# bshf: Non-relativistic B-spline Hartree-Fock

A program that solves the non-relativistic radial Schrödinger equation for atomic systems written in C++ for the first PHYS4070 assignment.

__Warning: Full functionality is limited to Lithium only! (for now...)__

## Dependencies:
* GNU GCC/G++ version >= 13.
* GNU Make version >= 3.80.
* LAPACK/BLAS version >= 3.10.
* libomp version >= 17.0.0.

## Compilation and usage
To build, run `make -j8` from the base directory.
The compiled program can then be run by calling `./bshf [args]`.
For documentation, use `./bshf --help`.
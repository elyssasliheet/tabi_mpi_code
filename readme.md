### This the repository contains Fortran code for the treecode-accelerated Poisson-Boltzmann Solver optimized with a preconditioner and MPI-based parallelization.

## by Elyssa Sliheet

This software package employs a well-conditioned boundary integral formulation for the electrostatic potential and its normal derivative on the molecular surface. The surface is triangulated and the integral equations are discretized by centroid collocation. The linear system is solved by GMRES iteration and the matrix-vector product is carried out by a Cartesian terraced which reduces the cost from O(N^2) to O(N*logN), where N is the number of faces in the triangulation. The TABI solver can be applied to compute the electrostatic potential on molecular surface and solvation energy.

This software includes added optimaztion via MPI-based parallelization, both sequentially and cyclically.

Four areas of parallelization-
1. RHS computation 
2. Treecode calculation 
3. Preconditioner
4. Solvation energy calculation

Original code can be found on sourceforge (https://sourceforge.net/p/tabipb/code/ci/master/tree/) by my PhD advisor, Dr. Weihua Geng.

Three molecular surface generators used here is called MSMS susrface (https://ccsb.scripps.edu/msms/downloads/)

A diagonal-block conditioner (see the reference below) is used to speed-up and stablize ill-conditioning caused by triangulation quality.


HOW TO RUN (on SMU's Maneframe III):

- module load intel
- module load mpi
- make (after navigating to appropriate file)
- sbatch test_scripts/test_64.job (to submit a job with 64 mpi tasks)
- run plotting.py (to visualize code of CPU times for treecode calculations for cyclic vs sequential implementations)


REFERENCES:

W. Geng and R. Krasny, A treecode-accelerated boundary integral Poisson-Boltzmann solver for continuum electrostatics of solvated biomolecules, J. Comput. Phys., 247, 62-87 (2013).

J. Chen and W. Geng, On preconditioning the treecode-accelerated boundary integral (TABI) Poisson-Boltzmann solver, J. Comput. Phys., 373, 750-762 (2018).

The C++ version of TABI is mained by Dr. Leighton Wilson on Github (https://github.com/Treecodes/TABI-PB)
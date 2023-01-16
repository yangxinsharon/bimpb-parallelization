# References
Weihua Geng, Ferosh Jacob. (2013). “A GPU-accelerated direct-sum boundary integral
Poisson–Boltzmann solver”. Computer Physics Communications, 184(6), 1490-1496.

# Acknowledgement
NSF Fund: DMS-1819193, DMS-2110896
Southern Methodist University Math Graduate Course: Daniel Reynolds - 6370 Parallel Computing

# bimpb-parallelization
You will find a serial version, and four parallelized versions including OpenMP, MPI, CUDA, and Kokkos of the software "bimpb" here. The bimpb uses "msms", a software which outputs the surface triangulization, and "gmres", a package which solves the linear system Ax = b with preconditioning.

gmres - translated by f2c and developed by Univ. of Tennessee and Oak Ridge National Laboratory in 1993.
test_proteins - path to use and save .pqr and .xyzr files
In readin.c files, find "fpath" to tune the location of test_proteins.
In main*.c files, find "fname" and "density" to tune these two parameters. The default selection is "1ajj" and "1".

## Serial, OpenMP, MPI, CUDA and Kokkos
common files: gl_constants.h gl_functions.h gl_variables.h pp_timer.c pp_timer.h readin.c msms
The serial, OpenMP and MPI share the common files in the sub-directory "serial_omp_mpi" with different suffix.

The serial version aims to compute the coulombic potential energy by solving Poisson-Boltzmann equation using Boundary Integral Method (bim-pb). 
All key computations are located in matvec.c file, including void matvecmul(), void comp_pot(), and void comp_source(). After computing the source term "b", the "gmres" package is used to compute the vector "x", and then the potential energy is computed in "comp_pot()".

The OpenMP (omp) version is developed based on the serial version. The main_omp.c is the same as main.c. In the matvec_omp.c file, void matvecmul(), void comp_pot(), and void comp_source() are parallelized through "for loops" so that a huge amount of tasks is shared across threads. 

The MPI version uses all processors. In the main_mpi.c file, it has to initiate using "MPI_Init" and finalize using "MPI_Finalize()". Only the first processor (root) calls the "readin()", and then it broadcasts to all other processors. In the matvec_mpi.c file, "matvecmul()", "comp_pot()", and "comp_source()" is parallelized by chunking the whole task size into each processor's interval. 

The CUDA version is developed by Jiahui Chen, who graduated in 2019 from SMU. The huge tasks are computed on GPU. 

The Kokkos uses C++ language, so this version mix compiles C and C++. It also computes on GPU. 

Examples:
Serial: 
$ ./bimpb.exe 
$ ./bimpb.exe 1ajj 1
OpenMP:
Login to HPC like ManeFrame II.
$ export OMP_NUM_THREADS=4
$ ./bimpb_omp.exe
MPI:
Login to HPC like ManeFrame II.
$ salloc -p standard-mem-s -N8 -n256 --x11=first
$ ./bimpb_mpi.exe
CUDA:
The CUDA version is developed by Jiahui Chen, who graduated in 2019 from SMU. 
$ srun -p v100x8 --gres=gpu:1 ./bimpb_cuda.exe
Kokkos:
$ srun -p development -c 4 --mem=16G --gres=gpu:volta:1 --pty $SHELL

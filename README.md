# References
Weihua Geng, Ferosh Jacob. (2013). “A GPU-accelerated direct-sum boundary integral
Poisson–Boltzmann solver”. Computer Physics Communications, 184(6), 1490-1496.

# Acknowledgement
NSF Fund: DMS-1819193, DMS-2110896
Southern Methodist University Math Graduate Course: Daniel Reynolds - 6370 Parallel Computing

# bimpb-parallelization
You will find a serial version, and four parallelized versions including OpenMP, MPI, CUDA, and Kokkos of the software "bimpb" here. The bimpb uses "msms", a software which outputs the surface triangulization, and "gmres", a package which solves the linear system Ax = b with preconditioning.

gmres - translated by f2c and developed by Univ. of Tennessee and Oak Ridge National Laboratory in 1993. \
test_proteins - path to use and save .pqr and .xyzr files \
In readin.c files, find "fpath" to tune the location of test_proteins.

## Serial, OpenMP, MPI, CUDA and Kokkos
common files: gl_constants.h gl_functions.h gl_variables.h pp_timer.c pp_timer.h readin.c msms \
The serial, OpenMP and MPI share the common files in the sub-directory "serial_omp_mpi" with different suffix.

The serial version aims to compute the coulombic potential energy by solving Poisson-Boltzmann equation using Boundary Integral Method (bim-pb). \
All key computations are located in matvec.c file, including void matvecmul(), void comp_pot(), and void comp_source(). After computing the source term "b", the "gmres" package is used to compute the vector "x", and then the potential energy is computed in "comp_pot()".

The OpenMP (omp) version is developed based on the serial version. The main_omp.c is the same as main.c. In the matvec_omp.c file, void matvecmul(), void comp_pot(), and void comp_source() are parallelized through "for loops" so that a huge amount of tasks is shared across threads. 

The MPI version uses all processors. In the main_mpi.c file, it has to initiate using "MPI_Init" and finalize using "MPI_Finalize()". Only the first processor (root) calls the "readin()", and then it broadcasts to all other processors. In the matvec_mpi.c file, "matvecmul()", "comp_pot()", and "comp_source()" is parallelized by chunking the whole task size into each processor's interval. 

The CUDA version is developed by Jiahui Chen, who graduated in 2019 from SMU. The huge tasks are computed on GPU. 

The Kokkos uses C++ language, so this version mix compiles C and C++. It also computes on GPU. 

Examples:
In main*.c files, find "fname" and "density" to tune these two parameters. The default selection is "1ajj" and "1".

MPI: \
Login to HPC like ManeFrame II. \
$ salloc -p standard-mem-s -N8 -n256 --x11=first \
$ module load gcc-9.2 hpcx \
$ srun -n 8 ./bimpb_mpi.exe \

Serial: \
$ ./bimpb.exe (1ajj) (1) \

OpenMP: \
$ export OMP_NUM_THREADS=4 \
$ ./bimpb_omp.exe \

CUDA:  \
Login to HPC like ManeFrame II.  \
$ module load nvhpc-22.2  \
$ srun -p v100x8 --gres=gpu:1 ./bimpb_cuda.exe  \

Kokkos:  \
Login to HPC like ManeFrame II.  \
$ srun -p development -c 4 --mem=16G --gres=gpu:volta:1 --pty $SHELL \
$ module load spack gcc-9.2 \
$ . /hpc/spack/share/spack/setup-env.sh \
$ spack load kokkos/qu45u5v \
$ cmake . \
$ make \
$ ./bimpb_kokkos.exe \



## SMU SuperPOD (must use VPN):
$ ssh username@slogin-01.superpod.smu.edu \

CUDA (SuperPOD):  can make and can run \
$ module load dev \
$ module load cuda-11.4.4-gcc-10.3.0-ctldo35 \
$ srun -G 1 ./bimpb_cuda.exe 1ajj 1 \


Kokkos (SuperPOD): \
srun -N 1 -G 1 -c 128 --mem=128G --time=12:00:00 --pty $SHELL \
Update on Oct 19, 2023: currently SuperPOD has 3 versions: \
        kokkos/3.7.01-jzhgq6o (not make) \
        kokkos/3.7.01-6zpfzzw (return segmentation fault) \
        kokkos/4.0.01-el36ysw \
Load the following packages in order: \
module load gcc/11.2.0 \
module load cuda/11.8.0-vbvgppx \
module load kokkos/4.0.01-el36ysw \

## To use Kokkos Kernels on SuperPOD at SMU:
Load the packages in order: \
module load gcc/11.2.0 \
module load cmake/3.26.3-utseokk \
module load cuda/11.8.0-vbvgppx \
module load kokkos/4.0.01-el36ysw { module load kokkos/4.1.00-n7m5qva (segmentation fault) } \
module load kokkos-kernels/4.1.00-xo4ovdm \
module load kokkos-nvcc-wrapper/4.1.00-hbetgcu \
cmake . -DKokkosKernels_SOURCE_DIR=$HOME/repos/kokkos-kernels \

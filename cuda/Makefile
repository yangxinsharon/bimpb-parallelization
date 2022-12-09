# makefile for screened coulombic potential computation with cuda
comp = nvcc
flag = -c -O3 -arch=sm_61 -use_fast_math
# for -arch=sm_xx check http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list
bimpb_cuda.exe: main_cuda.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec_cuda.o pp_timer.o
	$(comp) -o bimpb_cuda.exe *.o -lm
matvec_cuda.o: matvec_cuda.cu
	$(comp) $(flag) matvec_cuda.cu
pp_timer.o: ../pp_timer.c
	$(comp) $(flag) ../pp_timer.c
main_cuda.o: main_cuda.c
	$(comp) $(flag) main_cuda.c
readin.o: readin.c
	$(comp) $(flag) readin.c
d_sign.o: ../gmres/d_sign.c
	$(comp) $(flag) ../gmres/d_sign.c
daxpy.o: ../gmres/daxpy.c
	$(comp) $(flag) ../gmres/daxpy.c
dcopy.o: ../gmres/dcopy.c
	$(comp) $(flag) ../gmres/dcopy.c
ddot.o: ../gmres/ddot.c
	$(comp) $(flag) ../gmres/ddot.c
dgemv.o: ../gmres/dgemv.c
	$(comp) $(flag) ../gmres/dgemv.c
dnrm2.o: ../gmres/dnrm2.c
	$(comp) $(flag) ../gmres/dnrm2.c
drot.o: ../gmres/drot.c
	$(comp) $(flag) ../gmres/drot.c
drotg.o: ../gmres/drotg.c
	$(comp) $(flag) ../gmres/drotg.c
dscal.o: ../gmres/dscal.c
	$(comp) $(flag) ../gmres/dscal.c
dtrsv.o: ../gmres/dtrsv.c
	$(comp) $(flag) ../gmres/dtrsv.c
gmres.o: ../gmres/gmres.c
	$(comp) $(flag) ../gmres/gmres.c
clean: 
	\rm -rf *.o
realclean: clean
	\rm -rf *.exe ../test_proteins/*.face ../test_proteins/*.vert *.face *.vert *~

###################################################################
#  Makefile for screened coulombic potential computation 
#  with treecode
#  Xin (Sharon) Yang
#  SMU Mathematics
#  bimpb with Weihua Geng, Jiahui Chen
###################################################################

# compilers
# CC = clang # on MacOS
CC = gcc # on Linux, i.e., HPC-M2
# MPICC = mpicc

# flags
CFLAGS = -c -O3
OMPFLAGS = -fopenmp

# makefile targets
all: bimpb_tc.exe #bimpb_omp.exe bimpb_mpi.exe 

bimpb_tc.exe: main.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec.o  utilities.o xmalloc.o pp_timer.o
	$(CC) -p -o $@ $^ -lm
# bimpb_omp.exe: main.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec_omp.o pp_timer.o
# 	$(CC) $(OMPFLAGS) -p -o $@ $^ -lm
# bimpb_mpi.exe: main_mpi.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec_mpi.o pp_timer.o 
# 	$(MPICC) -o $@ $^ -lm
treecode.o: treepart/treecode.c
	$(CC) $(CFLAGS) treepart/treecode.c
utilities.o: treepart/utilities.c
	$(CC) $(CFLAGS) treepart/utilities.c
xmalloc.o: treepart/xmalloc.c
	$(CC) $(CFLAGS) treepart/xmalloc.c

pp_timer.o: pp_timer.c
	$(CC) $(CFLAGS) pp_timer.c
readin.o: readin.c
	$(CC) $(CFLAGS) readin.c
matvec.o: matvec.c
	$(CC) $(CFLAGS) matvec.c
main.o: main.c
	$(CC) $(CFLAGS) main.c

# matvec_omp.o: matvec_omp.c
# 	$(CC) $(CFLAGS) $(OMPFLAGS) matvec_omp.c
# main_omp.o: main_omp.c
# 	$(CC) $(CFLAGS) $(OMPFLAGS) main_omp.c

# matvec_mpi.o: matvec_mpi.c
# 	$(MPICC) $(CFLAGS) matvec_mpi.c
# main_mpi.o: main_mpi.c
# 	$(MPICC) $(CFLAGS) main_mpi.c

d_sign.o: gmres/d_sign.c
	$(CC) $(CFLAGS) gmres/d_sign.c
daxpy.o: gmres/daxpy.c
	$(CC) $(CFLAGS) gmres/daxpy.c
dcopy.o: gmres/dcopy.c
	$(CC) $(CFLAGS) gmres/dcopy.c
ddot.o: gmres/ddot.c
	$(CC) $(CFLAGS) gmres/ddot.c
dgemv.o: gmres/dgemv.c
	$(CC) $(CFLAGS) gmres/dgemv.c
dnrm2.o: gmres/dnrm2.c
	$(CC) $(CFLAGS) gmres/dnrm2.c
drot.o: gmres/drot.c
	$(CC) $(CFLAGS) gmres/drot.c
drotg.o: gmres/drotg.c
	$(CC) $(CFLAGS) gmres/drotg.c
dscal.o: gmres/dscal.c
	$(CC) $(CFLAGS) gmres/dscal.c
dtrsv.o: gmres/dtrsv.c
	$(CC) $(CFLAGS) gmres/dtrsv.c
gmres.o: gmres/gmres.c
	$(CC) $(CFLAGS) gmres/gmres.c
clean: 
	\rm -rf *.o *.out
realclean: clean
	\rm -rf *.exe ../test_proteins/*.face ../test_proteins/*.vert *.face *.vert *~

####### End of Makefile #######

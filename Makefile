FC = gfortran
FLAGS = -O3 
BLASDIR = $(HOME)/BLAS-3.6.0
LAPACKDIR = $(HOME)/lapack-3.6.0
FFTDIR = $(HOME)/fftw/.libs

eigen_par.x: eigen_par.o
	$(FC) $(FLAGS) -o eigen_par.x eigen_par.o $(LAPACKDIR)/liblapack.a $(BLASDIR)/blas_LINUX.a $(FFTDIR)/libfftw3.a 

eigen_par.o: eigen_par.F
	$(FC) $(FLAGS) -I$(HOME)/fftw/include -c -o eigen_par.o eigen_par.F

clean:
	rm -f eigen_par.x eigen_par.o 

all: eigen_par.x
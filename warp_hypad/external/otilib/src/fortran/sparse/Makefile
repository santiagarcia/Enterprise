
# Set gfortran or ifort to compile
FC = ifort # Default compiler

# Flags:
#   -O3: Optimization level 3 for best performance.
#   -fPIC: For shared library.
#   -fopenmp: For compatibility with parallel calls from Abaqus and others.
FCFLAGS = -O3 -fPIC -fopenmp
LCFLAGS = -I../core/ -L../core/ -lmaster_parameters -lreal_utils


SOURCES = sparse_oti.f90 

all: $(subst .f90,.mod, $(SOURCES))

sparse_oti.mod: sparse_oti.f90
	$(FC) -c $*.f90 $(FCFLAGS) $(LCFLAGS)
	ar rc lib$*.a $*.o
	ranlib lib$*.a

changegfortran:
	$(eval FC = gfortran)

gfortran: changegfortran $(subst .f90,.mod, $(SOURCES))

clean:
	-rm -f *.o *.mod *.a

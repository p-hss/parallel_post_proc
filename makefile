sourcedir = src/
objectdir = obj/

vpath %.f90 $(sourcedir)
vpath %.o $(objectdir) 
vpath %.mod $(objectdir) 

ifeq ($(shell hostname),gamling)
CFlaggs = -J$(objectdir) -O3 
LIB = -llapack -lblas -L/home/gamling/p.hess/lib -L/home/gamling/p.hess/Downloads/lapack-3.8.0
FC = mpif90 
endif

ifeq ($(shell hostname),cluster-i)
CFlaggs = -J$(objectdir) -O3 
LIB = -llapack -lblas -L/home/gamling/p.hess/lib -L/home/gamling/p.hess/Downloads/lapack-3.8.0
FC = mpif90 
endif


ifeq ($(shell hostname),cluster-a)
CFlaggs = -J$(objectdir) -O3 
LIB = -llapack -lblas -L/home/gamling/p.hess/lib -L/home/gamling/p.hess/Downloads/lapack-3.8.0
FC = mpif90 
endif

ifeq ($(findstring draco,$(shell hostname)),draco)
CFlaggs = -J$(objectdir) -O3 
LIB = ${MKL_HOME}/lib/intel64/libmkl_blas95_ilp64.a ${MKL_HOME}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
FC = mpif90 -f90=gfortran
endif

objects = main.o line_evolution.o diagnose.o input.o output.o

run: $(patsubst %.o, $(objectdir)%.o, global.o $(objects))
	$(FC) $(CFlaggs) -o $@ $^ $(LIB)

global.mod: $(sourcedir)global.f90  
	$(FC) $(CFlaggs) -c $(sourcedir)global.f90 

$(objectdir)global.o: 
	$(FC) $(CFlaggs) -c $(sourcedir)global.f90 -o $(objectdir)global.o

$(objectdir)main.o: $(objectdir)global.mod 
	$(FC) $(CFlaggs) -c $(sourcedir)main.f90 -o $(objectdir)main.o

$(objectdir)line_evolution.o: $(objectdir)global.mod
	$(FC) $(CFlaggs) -c $(sourcedir)line_evolution.f90 -o $(objectdir)line_evolution.o

$(objectdir)diagnose.o: $(objectdir)global.mod
	$(FC) $(CFlaggs) -c $(sourcedir)diagnose.f90 -o $(objectdir)diagnose.o

$(objectdir)input.o: $(objectdir)global.mod 
	$(FC) $(CFlaggs) -c $(sourcedir)input.f90 -o $(objectdir)input.o

$(objectdir)output.o: $(objectdir)global.mod
	$(FC) $(CFlaggs) -c $(sourcedir)output.f90 -o $(objectdir)output.o

.PHONY : clean

clean:
	rm -f run 
	rm -f $(objectdir)*
	rm -f src/*.mod
	rm -f *.mod

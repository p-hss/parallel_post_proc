sourcedir = src/
objectdir = obj/

vpath %.f90 $(sourcedir)
vpath %.o $(objectdir) 
vpath %.mod $(objectdir) 

CFlaggs = -J$(objectdir) -O3 

ifeq ($(shell hostname),gamling)
endif

ifeq ($(shell hostname),cluster-i)
endif


ifeq ($(shell hostname),cluster-a)
endif

ifeq ($(findstring draco,$(shell hostname)),draco)
endif

LIB = -llapack -lblas -L/home/gamling/p.hess/lib -L/home/gamling/p.hess/Downloads/lapack-3.8.0
FC = mpif90 

objects = main.o line_evolution.o diagnose.o input.o output.o

run: $(patsubst %.o,$(objectdir)%.o,global.o $(objects))
	$(FC) $(CFlaggs) -o $@ $^ $(LIB)

global.mod: $(sourcedir)global.f90  
	$(FC) $(CFlaggs) -c $(sourcedir)global.f90 

$(objectdir)global.o: 
	$(FC) $(CFlaggs) -c $(sourcedir)global.f90 -o $(objectdir)global.o

#%.o: $(objectdir)global.mod %.f90
#	$(FC) $(CFlaggs) -c $(sourcedir)$(patsubst %.o,%.f90, $(notdir $@)) -o $@

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

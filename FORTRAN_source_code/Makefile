#system variables
#fflags versions with "-Dopenmp" generate OpenMP parallel executables while
#versions without "-Dopenmp" generate serial executables

# intel fortran compiler
#f90 = ifort
#fflags = -implicitnone -O3 -march=native -Difort
#fflags = -implicitnone -O3 -march=native -qopenmp -Dopenmp -Difort

# gfortran compiler
f90 = gfortran 
#fflags = -fimplicit-none -O3 -march=native 
fflags = -fimplicit-none -O3 -fopenmp -march=native -Dopenmp 

# nvfortran
#f90 = nvfortran
#fflags = -O3 -march=native  
#fflags = -O3 -march=native -mp -Dopenmp 

# should not need to edit below this line
object_list = reconstruct.o data.o med_stats.o \
correlation.o grad_solve.o calc_ensemble.o setup.o
module_list = correlation_mod.mod ensemble_mod.mod med_stats_mod.mod \
data_mod.mod grad_solve_mod.mod setup_mod.mod
rm_lst = $(object_list) $(module_list) reconstruct 

all : reconstruct

install : reconstruct
	cp reconstruct ../wrappers_and_examples

reconstruct : $(object_list)
	$(f90) -o $@ $^ $(fflags)

reconstruct.o : reconstruct.f90 data.o med_stats.o correlation.o grad_solve.o calc_ensemble.o setup.o
	$(f90) -c $(fflags) $<

setup.o : setup.f90 data.o med_stats.o correlation.o grad_solve.o calc_ensemble.o 
	$(f90) -c $(fflags) $<

data.o : data.f90
	$(f90) -c $(fflags) $<

med_stats.o : med_stats.f90 
	$(f90) -c $(fflags) $<

correlation.o : correlation.F90 data.o
	$(f90) -c $(fflags) $<

grad_solve.o : grad_solve.f90 correlation.o data.o
	$(f90) -c $(fflags) $<

calc_ensemble.o : calc_ensemble.F90 grad_solve.o med_stats.o correlation.o
	$(f90) -c $(fflags) $<

.PHONY: clean
clean: 
	-@rm -f $(rm_lst)




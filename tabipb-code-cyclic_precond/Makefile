#F90= ifort
#flag= -fast -c  
#F90 = gfortran
#LIBS = -llapack -lblas -lm
#$LIBS = -lm
#F90 = mpif90
F90 = mpif90
flag= -O3 -c 
tabipb.exe: var_modules.o treecode3d_pb.o main.o readin.o dgmres_dep.o dgetrf_dep.o
	$(F90) -O3 -p -o tabipb.exe *.o 
var_modules.o:	var_modules.f90
	$(F90) $(flag)  var_modules.f90
main.o:		main.f90
	$(F90) $(flag)  main.f90
readin.o:	readin.f90
	$(F90) $(flag) readin.f90
treecode3d_pb.o:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
dgmres_dep.o:	dgmres_dep.f
	$(F90) $(flag) dgmres_dep.f
dgetrf_dep.o:	dgetrf_dep.f
	$(F90) $(flag) dgetrf_dep.f
molecule.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
comdata.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
bicg.mod:   	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode.mod:	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode3d_procedures.mod:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
clean: 
	rm *.o *.mod tabipb.exe surface_potential.dat

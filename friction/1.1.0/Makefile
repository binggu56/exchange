objects = qm.o ran2.o  gasdev.o  fitp.o
FC = gfortran
switch = -O3
EXEC = qm
#LIB = /usr/lib64/atlas/liblapack.so.3.0
LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

#cdat.o:mod.f
#	$(FC) -c mod.f
#cdat.mod: mod.f
#	$(FC) -c mod.f
fitp.o:fitp.f90
	$(FC) -c $(LIB) fitp.f90
%.o: %.f 
	$(FC) -c  $<

%.o: %.f90 
	$(FC) -c  $<

clean:
	rm *.o 
veryclean:
	rm *.o *.mod *.dat

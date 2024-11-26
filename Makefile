#compiler
FC=gfortran

#flags for the compiler
FFLAGS=-g -cpp -fdefault-real-8 

#files needed for compilation

#SRC=misc.f90 recon.f90 eos.f90 hllc.f90 test_hllc.f90
#SRC = io.f90 test_io.f90
SRC = sim_data.f90 misc.f90 read_par.f90 grid_init.f90 eos.f90 \
			applyBC.f90 sim_init.f90 sim_restart.f90 io.f90 recon.f90 hllc.f90 rk2.f90 main.f90

#object files have .o extension
OBJ=${SRC:.f90=.o}

#Loop (%) over *.f90 files to create *.o files
#$@ represents the target files and $< represents the prerequisite files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#executable to be created
EXEC=run.exe
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@

#clean all the compiled file and executable
clean:
	rm -f *.mod *.o  
cleanall:
	rm -f *.mod *.o *.exe
cleandata:	
	rm *.dat

CC=gcc
CFLAGS= 
F90=gfortran
F90FLAGS=-O3  
PARAFLAGS=-fopenmp

SOURCES= rdf_lib.f90 main.f90 


OBJ=$(addsuffix .o, $(basename $(SOURCES)))


.SUFFIXES :.c .f90

.f90.o:
	$(F90) $(F90FLAGS) $(PARAFLAGS) -c  $< 

.c.o:
	$(CC) $(CFLAGS) -c $<



rdf: $(OBJ) 
	@echo "Building new RDF ... "
	@echo "Current objects: $(OBJ)"
	$(F90) -o $@ $(F90FLAGS) $(PARAFLAGS) $(OBJ) 


clean:	
	@echo "Cleaning Directory ... "
	rm -f $(OBJ) 
 
veryclean: 
	@echo "Cleaning Directory AND Executables ..."
	rm -f $(OBJ) classical_md *.mod 

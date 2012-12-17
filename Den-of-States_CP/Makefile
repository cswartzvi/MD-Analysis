CC=gcc
CFLAGS= 
F90=ftn
F90FLAGS= -O3  
#PARAFLAGS=-fopenmp

SOURCES=  main.f90 


OBJ=$(addsuffix .o, $(basename $(SOURCES)))


.SUFFIXES :.c .f90

.f90.o:
	$(F90) $(F90FLAGS) $(PARAFLAGS) -c  $< 

.c.o:
	$(CC) $(CFLAGS) -c $<



cp-dos: $(OBJ)
	@echo "Building cp-dos ... "
	@echo "Current objects: $(OBJ)"
	$(F90) $(F90FLAGS) $(PARAFLAGS) $(OBJ) -o $@


clean:	
	@echo "Cleaning Directory ... "
	rm -f $(OBJ) 
 
veryclean: 
	@echo "Cleaning Directory AND Executables ..."
	rm -f $(OBJ) cp-dos *.mod 

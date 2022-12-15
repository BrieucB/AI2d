# définition des cibles particulières

.PHONY: clean

# désactivation des règles implicites

.SUFFIXES:

CC = gcc
LDLIBS = -lm -fopenmp
CFLAGS = -Ofast -fopenmp
CCFLAGS = #-g -Wall #-pg

OBJ = solve_PDEs_1d.o

solve_PDEs_1d : $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDLIBS) $(CCFLAGS) 

solve_PDEs_1d.o : solve_PDEs_1d.c
	$(CC) $(CFLAGS) -c $< -o $@ $(CCFLAGS)

clean :
	rm $(OBJ)

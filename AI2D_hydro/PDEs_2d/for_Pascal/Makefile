# définition des cibles particulières

.PHONY: clean

# désactivation des règles implicites

.SUFFIXES:

CC = gcc
LDLIBS = -lm -fopenmp
CFLAGS = -Ofast -fopenmp
CCFLAGS = #-g -Wall #-pg

OBJ = solve_PDEs_2d.o allocation.o exit_if.o erreurs.o

solve_PDEs_2d : $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDLIBS) $(CCFLAGS) 

solve_PDEs_2d.o : solve_PDEs_2d.c
	$(CC) $(CFLAGS) -c $< -o $@ $(CCFLAGS)

allocation.o : allocation.c allocation.h
	$(CC) $(CFLAGS) -c $< -o $@

exit_if.o : exit_if.c exit_if.h
	$(CC) $(CFLAGS) -c $< -o $@

erreurs.o : erreurs.c erreurs.h
	$(CC) $(CFLAGS) -c $< -o $@

clean :
	rm $(OBJ)

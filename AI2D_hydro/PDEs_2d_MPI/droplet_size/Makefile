# définition des cibles particulières

.PHONY: clean

# désactivation des règles implicites

.SUFFIXES:

CC = mpicc
LDLIBS = -lm
CFLAGS = -Ofast
CCFLAGS = -g -Wall #-pg

OBJ = solve_PDEs_2d.o allocation.o exit_if.o erreurs.o tools.o

solve_PDEs_2d : $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDLIBS) $(CCFLAGS) 

solve_PDEs_2d.o : solve_PDEs_2d.c
	$(CC) $(CFLAGS) -c $< -o $@ $(CCFLAGS)

tools.o : tools.c tools.h
	$(CC) $(CFLAGS) -c $< -o $@

allocation.o : allocation.c allocation.h
	$(CC) $(CFLAGS) -c $< -o $@

exit_if.o : exit_if.c exit_if.h
	$(CC) $(CFLAGS) -c $< -o $@

erreurs.o : erreurs.c erreurs.h
	$(CC) $(CFLAGS) -c $< -o $@

clean :
	rm $(OBJ)

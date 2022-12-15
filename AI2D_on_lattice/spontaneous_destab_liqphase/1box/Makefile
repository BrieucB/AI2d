# définition des cibles particulières

.PHONY: clean

# désactivation des règles implicites

.SUFFIXES:

CC = gcc
LDLIBS = -lm #-fopenmp
CFLAGS = -Ofast #-fopenmp
CCFLAGS = #-g -Wall #-pg

OBJ = tools.o activeIsing.o allocation.o exit_if.o erreurs.o pcg_basic.o

activeIsing : $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) -o $@ $(LDLIBS) $(CCFLAGS) 

activeIsing.o : activeIsing.c
	$(CC) $(CFLAGS) -c $< -o $@ $(CCFLAGS)

tools.o : tools.c tools.h
	$(CC) $(CFLAGS) -c $< -o $@ $(CCFLAGS)

allocation.o : allocation.c allocation.h
	$(CC) $(CFLAGS) -c $< -o $@

exit_if.o : exit_if.c exit_if.h
	$(CC) $(CFLAGS) -c $< -o $@

erreurs.o : erreurs.c erreurs.h
	$(CC) $(CFLAGS) -c $< -o $@

pcg_basic.o : pcg_basic.c pcg_basic.h
	$(CC) $(CFLAGS) -c $< -o $@ 

clean :
	rm $(OBJ)

CC=gcc
COMPILER_FLAGS = -Wall -g -O3 -march=native -mpc64  -fno-exceptions -ffast-math -funroll-loops -fstrict-aliasing -finline-functions
CFLAGS=-I. $(COMPILER_FLAGS)

DEPS = errors.h io.h tools.h mol.h msd.h linreg.h main.h
OBJ = errors.o io.o tools.o mol.o msd.o linreg.o main.o

LIBS=-lm

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

anal: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *~ core $(INCDIR)/*~

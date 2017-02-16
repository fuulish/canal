CC = gcc
COMPILER_FLAGS = -Wall -g -O3 -march=native -mpc64  -fno-exceptions -ffast-math -funroll-loops -fstrict-aliasing -finline-functions
#COMPILER_FLAGS = -w -g -O3 -ffast-math
#COMPILER_FLAGS = -w -g -DOPENMP -fopenmp
LIBS = -lm
LINKER_FLAGS = $(LIBS)
OBJ_NAME = anal

OBJS = errors.c io.c tools.c mol.c msd.c main.c

all: $(OBJS) 
	$(CC) $(OBJS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)

clean:
	- rm $(OBJ_NAME)

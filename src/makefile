CC = g++
COMPILER_FLAGS = -w -g
LINKER_FLAGS = 
OBJ_NAME = anal

OBJS = errors.c io.c tools.c main.c

all: $(OBJS) 
	$(CC) $(OBJS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)

clean:
	- rm $(OBJ_NAME)

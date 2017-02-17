IDIR =../include
CC = gcc
COMPILER_FLAGS = -Wall -g -O3 -march=native -mpc64  -fno-exceptions -ffast-math -funroll-loops -fstrict-aliasing -finline-functions
CFLAGS = -I$(IDIR) $(COMPILER_FLAGS)

ODIR = obj
LDIR =../lib

LIBS = -lm

_DEPS = errors.h io.h tools.h mol.h msd.h linreg.h data.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = errors.o io.o tools.o mol.o msd.o linreg.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(ODIR) anal

anal: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR):
	mkdir $(ODIR)

.PHONY: clean distclean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

distclean: clean
	rm -rf $(ODIR) anal

tests: normal split

normal:
	./anal
	diff cond_all.out ../refs/cond_all.dat

split:
	./anal split
	diff cond_neinst.out ../refs/cond_neinst.dat
	diff cond_catcat.out ../refs/cond_catcat.dat
	diff cond_anicat.out ../refs/cond_anicat.dat
	diff cond_aniani.out ../refs/cond_aniani.dat

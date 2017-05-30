IDIR =../include
CC = gcc
COMPILER_FLAGS = -Wall -g -O3 -march=native -mpc64  -fno-exceptions -ffast-math -funroll-loops -fstrict-aliasing -finline-functions
#COMPILER_FLAGS = -Wall -g -O2 -march=native -mpc64  -fno-exceptions -funroll-loops -fstrict-aliasing -finline-functions
#CC = clang
#COMPILER_FLAGS = -Wall -g -O3 -march=native -fno-exceptions -ffast-math -funroll-loops -fstrict-aliasing -finline-functions
DFLAGS = -DDEBUG
CFLAGS = -I$(IDIR) $(COMPILER_FLAGS) $(DFLAGS)

ODIR = obj
LDIR =../lib
DDIR = ../data

LIBS = -lm -lgsl -lgslcblas

_DEPS = errors.h io.h tools.h mol.h vel.h msd.h linreg.h constants.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = errors.o io.o tools.o mol.o vel.o msd.o linreg.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(ODIR) canal

canal: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR):
	mkdir $(ODIR)

.PHONY: clean distclean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

distclean: clean
	rm -rf $(ODIR) canal *.out

#_INP = normal split spatial
#INP = $(patsubst %,$(DDIR)/%, $(_INP))
#
#tests:
#	./canal $@

tests: normal.out split.out spatial.out

normal.out: canal
	./canal ../data/normal.inp &> normal.out
	diff cond_all.out ../refs/cond_all.dat

split.out: canal
	./canal ../data/split.inp &> split.out
	diff cond_neinst.out ../refs/cond_neinst.dat
	diff cond_catcat.out ../refs/cond_catcat.dat
	diff cond_anicat.out ../refs/cond_anicat.dat
	diff cond_aniani.out ../refs/cond_aniani.dat

spatial.out: canal
	./canal ../data/spatial.inp &> spatial.out
	diff cond_neinst.out ../refs/cond_neinst_spatial.dat
	diff cond_catcat.out ../refs/cond_catcat_spatial.dat
	diff cond_anicat.out ../refs/cond_anicat_spatial.dat
	diff cond_aniani.out ../refs/cond_aniani_spatial.dat

check:	tests
	../tools/conform.sh

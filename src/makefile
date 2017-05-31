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

.PHONY: clean distclean check check-normal check-split check-spatial check-elmo

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

normal.out: all
	./canal ../data/normal.inp > normal.out

check-normal: normal.out
	diff cond_all.out ../refs/cond_all.dat

split.out: all
	./canal ../data/split.inp > split.out

check-split: split.out
	diff cond_neinst.out ../refs/cond_neinst.dat
	diff cond_catcat.out ../refs/cond_catcat.dat
	diff cond_anicat.out ../refs/cond_anicat.dat
	diff cond_aniani.out ../refs/cond_aniani.dat

spatial.out: all
	./canal ../data/spatial.inp > spatial.out
	mv cond_neinst.out cond_neinst_spatial.out
	mv cond_catcat.out cond_catcat_spatial.out
	mv cond_anicat.out cond_anicat_spatial.out
	mv cond_aniani.out cond_aniani_spatial.out

check-spatial: spatial.out
	diff cond_neinst_spatial.out ../refs/cond_neinst_spatial.dat
	diff cond_catcat_spatial.out ../refs/cond_catcat_spatial.dat
	diff cond_anicat_spatial.out ../refs/cond_anicat_spatial.dat
	diff cond_aniani_spatial.out ../refs/cond_aniani_spatial.dat

elmo.out: all
	./canal ../data/elmo.inp > elmo.out

check-elmo: elmo.out

check: check-normal check-split check-spatial check-elmo
	../tools/conform.sh

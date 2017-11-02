IDIR =../include
CC = gcc
COMPILER_FLAGS = -Wall -g -O3 -march=native -mpc64  -fno-exceptions -ffast-math -funroll-loops -fstrict-aliasing -finline-functions -std=c99
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

all: canal

canal: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(OBJ):	| $(ODIR)

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

normal.out: canal ../data/normal.inp
	./canal ../data/normal.inp > normal.out

check-normal: normal.out
	diff cond_all.out ../refs/cond_all.dat

split.out: canal ../data/split.inp
	./canal ../data/split.inp > split.out
	mv cond_neinaa.out cond_neinaa_split.out
	mv cond_neincc.out cond_neincc_split.out
	mv cond_catcat.out cond_catcat_split.out
	mv cond_anicat.out cond_anicat_split.out
	mv cond_aniani.out cond_aniani_split.out

check-split: split.out
	diff cond_neinaa_split.out ../refs/cond_neinaa_split.dat
	diff cond_neincc_split.out ../refs/cond_neincc_split.dat
	diff cond_catcat_split.out ../refs/cond_catcat_split.dat
	diff cond_anicat_split.out ../refs/cond_anicat_split.dat
	diff cond_aniani_split.out ../refs/cond_aniani_split.dat

spatial.out: canal ../data/spatial.inp
	./canal ../data/spatial.inp > spatial.out
	mv cond_neinaa.out cond_neinaa_spatial.out
	mv cond_neincc.out cond_neincc_spatial.out
	mv cond_catcat.out cond_catcat_spatial.out
	mv cond_anicat.out cond_anicat_spatial.out
	mv cond_aniani.out cond_aniani_spatial.out

check-spatial: spatial.out
	diff cond_neinaa_spatial.out ../refs/cond_neinaa_spatial.dat
	diff cond_neincc_spatial.out ../refs/cond_neincc_spatial.dat
	diff cond_catcat_spatial.out ../refs/cond_catcat_spatial.dat
	diff cond_anicat_spatial.out ../refs/cond_anicat_spatial.dat
	diff cond_aniani_spatial.out ../refs/cond_aniani_spatial.dat

elmo.out: canal ../data/elmo.inp
	./canal ../data/elmo.inp > elmo.out

check-elmo: elmo.out

check: check-normal check-split check-spatial check-elmo
	../tools/conform.sh

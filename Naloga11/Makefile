CC=clang
SRCNAM=cluster.c
PRGNAM=ccluster
OPTFLAGS=-march=corei7 -funroll-loops -ffast-math
WFLAGS=-Wall

all:
	$(CC) $(SRCNAM) -o $(PRGNAM) $(WFLAGS) $(OPTFLAGS)

clean:
	rm -f $(PRGNAM)

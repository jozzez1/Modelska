CC=gcc47
SRCNAM=cluster.c
PRGNAM=ccluster
OPTFLAGS=-march=corei7 -funroll-loops -ffast-math -O3
WFLAGS=-Wall

all:
	$(CC) $(SRCNAM) -o $(PRGNAM) $(WFLAGS) $(OPTFLAGS)

random:
	gcc47 randomize.c -o crandomize $(WFLAGS) $(OPTFLAGS)

clean:
	rm -f $(PRGNAM) crandomize

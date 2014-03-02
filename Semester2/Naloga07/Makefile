CC=		clang
PRGNAM=		prg
SRCNAM=		main.c
OPTFLAG=	-O3
CPUOPT=		-march=corei7 -ffast-math
LIBFLAGS=	-L/usr/local/lib -lm
INCLUDE=	-I/usr/local/include
WARN=		-Wall -pedantic

all:
	$(CC) -o $(PRGNAM) $(SRCNAM) $(OPTFLAG) $(CPUOPT) $(INCLUDE) $(LIBFLAGS) $(WARN) -g -pg

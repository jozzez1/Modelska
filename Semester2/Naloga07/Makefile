CC=		clang
PRGNAM=		prg
SRCNAM=		main.c
OPTFLAG=	-O3
CPUOPT=		-march=corei7
LIBFLAGS=	-L/usr/local/lib -lm -lgsl -lgslcblas
INCLUDE=	-I/usr/local/include
WARN=		-Wall

all:
	$(CC) -o $(PRGNAM) $(SRCNAM) $(OPTFLAG) $(CPUOPT) $(INCLUDE) $(LIBFLAGS) $(WARN) -g

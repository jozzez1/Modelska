CC=			clang
CPUFLAGS=	-march=corei7
OFLAG=		-O3
OPTFLAGS=	${OFLAG} ${CPUFLAG} \
			-ffast-math -fno-stack-protector \
			-funroll-loops
WARNFLAGS=	-Wall
LIBFLAGS=	-L/usr/local/lib \
			-lm -lgsl -lgslcblas
CFLAGS=		-I/usr/local/include \
			${OPTFLAGS} ${WARNFLAGS} -g
PRGNAM=		prg

all: main

main:
	${CC} main.c -o ${PRGNAM} ${CFLAGS} ${LIBFLAGS}

remove:
	rm -f ${PRGNAM}

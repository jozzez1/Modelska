CC=			clang
CPUFLAGS=	-march=corei7 -msse4.1 -mssse3
OFLAG=		-Ofast
OPTFLAGS=	${OFLAG} ${CPUFLAGS} \
			-ffast-math -fno-stack-protector \
			-funroll-loops
WARNFLAGS=	-Wall -Wno-typedef-redefinition
LIBFLAGS=	-L/usr/local/lib \
			-lm -lgsl -lgslcblas \
			-lmgl -lmgl-wnd
CFLAGS=		-I/usr/local/include \
			${OPTFLAGS} ${WARNFLAGS} -g
PRGNAM=		prg

all: main

main:
	${CC} main.c -o ${PRGNAM} ${CFLAGS} ${LIBFLAGS}

clean:
	rm -f anim/*

remove:
	rm -f ${PRGNAM}

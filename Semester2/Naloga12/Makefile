CC=			clang
PRGNAM=		prg
CPUOPT=		-march=corei7 -mssse3 -msse4.1
OFLAG=		-Ofast
OPTFLAGS=	${OFLAG} ${CPUOPT} \
			-funroll-loops -ffast-math -fno-stack-protector
WARN=		-Wall -Wno-return-type-c-linkage
LIBFLAGS=	-L/usr/local/lib \
	       	-lm
CFLAGS=		-I/usr/local/include \
			${WARN} ${OPTFLAGS} -g

all: main.o
	${CC} ${CFLAGS} ${LIBFLAGS} main.o navier_stokes.o sor.o -o ${PRGNAM}

main.o: navier_stokes.o sor.o
	${CC} ${CFLAGS} -c main.c

navier_stokes.o: sor.o
	${CC} ${CFLAGS} -c navier_stokes.c

sor.o:
	${CC} ${CFLAGS} -c sor.c

clean:
.for file in \
	sor.o \
	navier_stokes.o \
	prg.core
	rm -f ${file}
.endfor

remove:
	rm -f ${PRGNAM}

CC=		clang
PRGNAM=		prg
CPUOPT=		-march=corei7 -mssse3 -msse4.1
OFLAG=		-Ofast
OPTFLAGS=	${OFLAG} ${CPUOPT} \
		-funroll-loops -ffast-math -fno-stack-protector
WARN=		-Wall
LIBFLAGS=	-L/usr/local/lib \
	       	-lm
CFLAGS=		-I/usr/local/include \
		${WARN} ${OPTFLAGS} -g

all: main.o global.o
	${CC} ${CFLAGS} ${LIBFLAGS} main.o navier_stokes.o sor.o global.o -o ${PRGNAM}

main.o: navier_stokes.o global.o
	${CC} ${CFLAGS} -c main.c

navier_stokes.o: sor.o global.o
	${CC} ${CFLAGS} -c navier_stokes.c

sor.o: global.o
	${CC} ${CFLAGS} -c sor.c

global.o:
	${CC} ${CFLAGS} -c global.c

clean:
.for file in \
	main.o \
	navier_stokes.o \
	sor.o \
	global.o \
	prg.gmon \
	prg.core
	rm -f ${file}
.endfor

remove:
	rm -f ${PRGNAM}

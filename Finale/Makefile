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
		${WARN} ${OPTFLAGS} -pipe -g -pg

all: main.o
	${CC} ${CFLAGS} ${LIBFLAGS} main.o zero_solve.o star_system.o -o ${PRGNAM}

main.o: zero_solve.o star_system.o
	${CC} ${CFLAGS} -c main.c

star_system.o: zero_solve.o
	${CC} ${CFLAGS} -c star_system.c

zero_solve.o:
	${CC} ${CFLAGS} -c zero_solve.c

purge: clean remove
	rm -f prg.gmon

clean:
.for file in \
	main.o \
	star_system.o \
	zero_solve.o \
	prg.core
	rm -f ${file}
.endfor

remove:
	rm -f ${PRGNAM}

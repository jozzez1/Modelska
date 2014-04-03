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
		${WARN} ${OPTFLAGS} -g -pg

all: main.o
	${CC} ${CFLAGS} ${LIBFLAGS} main.o zero_solve.o -o ${PRGNAM}

main.o: zero_solve.o
	${CC} ${CFLAGS} -c main.c

zero_solve.o:
	${CC} ${CFLAGS} -c zero_solve.c

purge: clean remove
	rm -f prg.gmon

clean:
.for file in \
	main.o \
	zero_solve.o \
	prg.core
	rm -f ${file}
.endfor

remove:
	rm -f ${PRGNAM}

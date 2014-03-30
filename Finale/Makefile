CXX=		clang++
PRGNAM=		prg
CPUOPT=		-march=corei7 -mssse3 -msse4.1
OFLAG=		-Ofast
OPTFLAGS=	${OFLAG} ${CPUOPT} \
		-funroll-loops -ffast-math -fno-stack-protector
WARN=		-Wall -Wno-return-type-c-linkage
LIBFLAGS=	-L/usr/local/lib \
	       	-lm
CXXFLAGS=	-I/usr/local/include \
		${WARN} ${OPTFLAGS}

all: main.o global.o
	${CXX} ${CXXFLAGS} ${LIBFLAGS} main.o global.o trotter.o -o ${PRGNAM}

main.o: trotter.o global.o
	${CXX} ${CXXFLAGS} -c main.cxx

trotter.o: global.o
	${CXX} ${CXXFLAGS} -c trotter.cxx

global.o:
	${CXX} ${CXXFLAGS} -c global.cxx

clean:
.for file in \
	main.o \
	trotter.o \
	global.o \
	prg.core
	rm -f ${file}
.endfor

remove:
	rm -f ${PRGNAM}

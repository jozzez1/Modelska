CXX=		clang++
PRGNAM=		prg
CPUOPT=		-march=corei7 -mssse3 -msse4.1
OFLAG=		-O3
OPTFLAGS=	${OFLAG} ${CPUOPT} \
		-funroll-loops -ffast-math -fno-stack-protector
WARN=		-Wall -Wno-return-type-c-linkage
LIBFLAGS=	-L/usr/local/lib \
	       	-lm -lmgl -lmgl-wnd
CXXFLAGS=	-I/usr/local/include/eigen3 \
		-I/usr/local/include \
		${WARN} ${OPTFLAGS} -g

all: main.o
	${CXX} ${CXXFLAGS} ${LIBFLAGS} main.o walker.o profiles.o -o ${PRGNAM}

main.o: walker.o
	${CXX} ${CXXFLAGS} -c main.cxx

walker.o: profiles.o
	${CXX} ${CXXFLAGS} -c walker.cxx

profiles.o:
	${CXX} ${CXXFLAGS} -c profiles.cxx

clean:
.for file in \
	main.o \
	walker.o \
	profiles.o \
	prg.core
	rm -f ${file}
.endfor

remove:
	rm -f ${PRGNAM}

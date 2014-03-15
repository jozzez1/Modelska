CXX=		clang++
PRGNAM=		prg
CPUOPT=		-march=corei7
OFLAG=		-O3
OPTFLAGS=	${OFLAG} ${CPUOPT} -ffast-math -funroll-loops
WARN=		-Wall
LIBFLAGS=	-L/usr/local -lm
CXXFLAGS=	-I/usr/local/include/eigen3 \
		${WARN} ${OPTFLAGS}

all: main.o walker.o
	${CXX} ${CXXFLAGS} ${LIBFLAGS} main.o walker.o profiles.o -o ${PRGNAM}

main.o: walker.o
	${CXX} ${CXXFLAGS} -c main.cxx

walker.o: profiles.o
	${CXX} ${CXXFLAGS} -c walker.cxx

profiles.o:
	${CXX} ${CXXFLAGS} -c profiles.cxx

clean:
	@rm -f main.o walker.o profiles.o

remove:
	rm -f ${PRGNAM}

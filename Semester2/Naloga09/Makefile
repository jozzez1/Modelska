CXX=clang++
PRGNAM=		prg
CPUOPT=	-march=corei7
OFLAG=		-O3
OPTFLAGS=	${OFLAG} ${CPUOPT} -ffast-math -funroll-loops
WARN=		-Wall
LIBFLAGS=	-L/usr/local -lm
CXXFLAGS=	-I/usr/local/include/eigen3 \
		${WARN} ${OPTFLAGS}

all: main.o walker.o
	$(CXX) -o $(PRGNAM) main.o walker.o $(CXXFLAGS) $(LIBFLAGS)

main.o:
	${CXX} ${CXXFLAGS} -c main.cxx

walker.o:
	${CXX} ${CXXFLAGS} -c walker.cxx

clean:
	rm ${PRGNAM} main.o walker.o

CXX=		clang++
CPUFLAG=	-march=corei7
OFLAG=		-O3
OPTFLAGS=	${OFLAG} ${CPUFLAG} \
			-ffast-math -fno-stack-protector \
			-funroll-loops
WARNFLAGS=	-Wall -Wno-return-type-c-linkage
LIBFLAGS=	-L/usr/local/lib \
			-lmgl -lmgl-wnd
CXXFLAGS=	-I/usr/local/include \
			${OPTFLAGS} ${WARNFLAGS} -g
PRGNAM=		prg

all: main

main:
	${CXX} main.cxx -o ${PRGNAM} ${CXXFLAGS} ${LIBFLAGS}

remove:
	rm -f ${PRGNAM}

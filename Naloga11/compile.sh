#!/bin/tcsh

set libdir=/usr/local/lib
set incdir=/usr/local/include
set prgnam=prg
set cflags=(-Wall -O2 -lgsl -lgslcblas -std=c99)

gcc -I $incdir -L $libdir main.c -o $prgnam $cflags

exit 0


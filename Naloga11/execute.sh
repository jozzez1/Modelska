#!/bin/tcsh

set dim=($1)
set con=($2)
set imax=($3)
set seed=0

if ( $# == 4 ) then
	seed=$4
endif

./prg $dim $con $imax $seed

set mfile="matrix-$1-$2.txt"
set ofile="connectivity-$1-$2.txt"

exit 0

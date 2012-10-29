#!/bin/tcsh

set dim=($1)
set con=($2)
set seed=0

if ( $# == 3 ) then
	set seed=($3)
endif

./prg $dim $con $seed

set mfile="matrix-$1-$2.txt"
set ofile="connectivity-$1-$2.txt"

exit 0

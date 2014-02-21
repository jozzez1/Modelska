#!/bin/sh

input=`echo $1 | sed 's/.pgm//g'`

convert $input.pgm test.txt
./reshaper test.txt > $input.dat

rm test.txt

# here comes the octave script
octave -q << INPUT

M=$2;

% load the data samples
u1 = load('-ascii', '$input.dat');

[Y, X] = size (u1);
[U, S, V] = svd (u1, 1);
V = V';

U = resize (U, Y, M);
S = resize (S, M, M);
V = resize (V, M, X);

u1 = U * S * V;

% now we have to make sure everything is rescaled
% in [0, 255] interval of gray ...
u1 = round (u1);
u1 = min(u1, 255);
u1 = max(u1, 0);

fid = fopen ('small-$input-M$2.pgm', 'w');
fprintf (fid, "P2\n%d %d\n255\n", X, Y);
for i = 1:Y
	fprintf (fid, "%d   ", u1(i, :));
	fprintf (fid, "\n");
end
fclose(fid);

INPUT

convert small-$input-M$2.pgm small-$input-M$2.png
viewnior small-$input-M$2.png &
rm small-$input-M$2.pgm

exit 0


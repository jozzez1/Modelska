#!/bin/sh

input=`echo $1 | sed 's/.pgm//g'`

convert $input.pgm test.txt
./reshaper test.txt > $input.dat

rm test.txt

# here comes the octave script
octave -q << INPUT

filter=$2;

% load the data samples
L0 = load('-ascii', '$input.dat');

% we have to arrange it to 1D
L0 = reshape (L0, 313*256, 1);
L0 = L0;

% FFT
fL0 = fft(L0);

% now it's time for deconvolution
t = (0:(313*256 - 1))';
v = exp(-t./30)./30;
fv = fft(v);

if filter == 1
	% Wiener filter
	Lc = load('-ascii', 'lincoln_L30_N00.dat');
	Lc = reshape (Lc, 313*256, 1);
	fLc = fft (Lc);
	n = fL0 - fLc;
	f = fLc;
	W = 1 ./ (1 .+ (abs(n) ./ abs(f)).^2);

	%let's add a bit of hamming
	H = 1 - hamming(313*256);
	W = H .* W;
else
	W = 1;
endif

u1 = real(ifft(W .* fL0 ./ fv));
u1 = reshape (u1, 256, 313);

% now we have to make sure everything is rescaled
% in [0, 255] interval of gray ...
u1 = round (u1);
u1 = min(u1, 255);
u1 = max(u1, 0);

u1 = reshape (u1', 16, 5008)';
fid = fopen ('fixed-$input-filter$2.pgm', 'w');
fprintf (fid, "P2\n313 256\n255\n");
for i = 1:5008
	fprintf (fid, "%d   ", u1(i, :));
	fprintf (fid, "\n");
end
fclose(fid);

INPUT

convert fixed-$input-filter$2.pgm fixed-$input-filter$2.png
viewnior fixed-$input-filter$2.png &
rm fixed-$input-filter$2.pgm

exit 0

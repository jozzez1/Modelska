% load the data samples
L0 = load('-ascii', 'lincoln_L30_N30.dat');

% we have to arrange it to 1D
L0 = reshape (L0, 313*256, 1);

% FFT
fL0 = fft(L0);

% now it's time for deconvolution
t = (0:(313*256 - 1))';
%t = flipud (t);
v = exp(-t./30)./30;
fv = fft(v);

% Wiener filter
Lc = load('-ascii', 'lincoln_L30_N30.dat');
Lc = reshape (Lc, 313*256, 1);
fLc = fft (Lc);
n = fL0 - fLc;
f = fLc;
W = 1 ./ (1 .+ (abs(n) ./ abs(f)).^2);

u1 = real(ifft(W .* fL0 ./ fv));
u1 = reshape (u1, 256, 313);

% now we have to make sure everything is rescaled
% in [0, 255] interval of gray ...
u1 = round (u1);
u1 = min(u1, 255);
u1 = max(u1, 0);

u1 = reshape (u1', 16, 5008)';
fid = fopen ('test1.pgm', 'w');
fprintf (fid, "P2\n313 256\n255\n");
for i = 1:5008
	fprintf (fid, "%d   ", u1(i, :));
	fprintf (fid, "\n");
end
fclose(fid);

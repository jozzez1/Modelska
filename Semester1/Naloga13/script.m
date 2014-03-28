%function [u0, u1, u2, u3, t] = compute ()
	t = 0:511;
	t = t.';
	v = (0.5/16) * exp (-(abs(t - 255)./16)); % prenosna funkcija

	% naše vhodne datoteke
	c0 = load('-ascii', 'signal0.dat');
	c1 = load('-ascii', 'signal1.dat');
	c2 = load('-ascii', 'signal2.dat');
	c3 = load('-ascii', 'signal3.dat');

	% izračunajmo transformiranke
	C0 = fft(c0);
	C1 = fft(c1);
	C2 = fft(c2);
	C3 = fft(c3);

	V  = fft(v);

	% izračunamo šum
	N1 = C0 - C1;
	N2 = C0 - C2;
	N3 = C0 - C3;

	% izračunamo signal
	F0 = C0;
	F1 = C1 - N1;
	F2 = C2 - N2;
	F3 = C3 - N3;

	% dodamo Wienerjev filter
	W0 = ones(512,1);
	W1 = 1 ./ (1 .+ (abs(N1)./abs(F1)).^2);
	W2 = 1 ./ (1 .+ (abs(N2).^2)./(abs(F2).^2));
	W3 = 1 ./ (1 .+ (abs(N3).^2)./(abs(F3).^2));

	% signal obrnemo nazaj
	u0 = ifft(W0 .* (C0 ./ V));
	u1 = ifft(W1 .* (C1 ./ V));
	u2 = ifft(W2 .* (C2 ./ V));
	u3 = ifft(W3 .* (C3 ./ V));

%	return;
%endfunction

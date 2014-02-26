function [P, z, A, sigma2] = spect (filename, order)
	% supress some annoying warnings
	warning ("off", "Octave:broadcast");
	warning ("off", "Octave:shadowed-function");
	pkg load tsa
	
	% load the data
	y = load ('-ascii', filename);
	
	% some spectrum parameters
	N = length(y);
	x = 0:N-1;
	w = exp(2*pi*i*x/N);
	
	% get the coefficients
	[A, RC, PE] = lattice (y.', order);
	
	% A = [Ap, ... A2, A1], the polynomial is 
	% p(x;A) = (-A_p) + (-A_{p-1}) x + (-A_{p-2}) x^2 + ... + (-A_1) x^{p-1} + x^p
	% so we have to fix this. We can use `ar2poly' function, or we can just
	% do this, for the sake of transparency:
	A = [1, -A];
	
	% we look for zeros of this polynomial
	z = roots (A);
	
	% now fix the divergent roots -- just in case
	N = length(z);
	k = 0;
	for i = 1:N
		if abs(z(i)) > 1
			z(i) = z(i)/abs(z(i));
			k++;
		endif
	end
	k
	
	% polonom values
	N = length(y);	
	p = -(z .- w);
	p = prod (p);
	
	sigma2 = max (PE);
	
	% and now the spectrum
	P = sigma2 ./ (abs(p).^2);

	% add those warnings back
	warning ("on", "Octave:broadcast");
	warning ("on", "Octave:shadowed-function");
	pkg unload tsa

	return;
endfunction

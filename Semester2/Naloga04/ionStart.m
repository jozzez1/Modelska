function [E, u, U] = ionStart (rmin, rmax, N)
	Z = 2; % (Helium atom)
	r = linspace (rmin, rmax, N);

	% this function is going to give us the starting point
	% for the iteration (much like the DFT)
	Zstar = Z - 5/16;

	u = 2*sqrt(Zstar) * Zstar * r .* exp (-Zstar * r);
	u = (u ./ sqrt(sumsq(u))) .* (sqrt(N / (rmax - rmin)));
	U = -potentialU (u, rmin, rmax, N);
	E = -1;

	return;
endfunction

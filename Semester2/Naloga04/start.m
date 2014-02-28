function [E, u, U] = start (E0, D, ell, rmin, rmax, N)
	r = linspace (rmin, rmax, N);
	V = -1 ./ r + (ell * (ell + 1)) ./ (2 * r .^ 2);

	if rmin == 0
		V (1) = 0;
	endif

	[E, u] = ground (V, E0, D, rmin, rmax, N);
	U = -potentialU (u, rmin, rmax, N);
	E *= 2;

	return;
endfunction

function [E, u, U] = step (R, phi, E0, D, rmin, rmax, N)
	% R -> u
	% V -> U

	% potential
	r = linspace (rmin, rmax, N);
	V = -2 ./ r - phi ./ r;

	if rmin == 0
		V(1) = 0;
	endif

	% first we need the new R, which is u
	[E, u] = ground (V, E0, D, rmin, rmax, N);

	% and now we calculate the new phi
	U = -potentialU (u, rmin, rmax, N);

	return;
endfunction

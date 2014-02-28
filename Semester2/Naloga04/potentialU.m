function U = potentialU (u, rmin, rmax, N)

	h = (rmax - rmin) / N;

	r = linspace (rmin, rmax, N);
	ur = - (u.^2) ./ r;
	ur(1) = 0;
	L = mLaplace (rmin, rmax, N);
	
	[nr, nc] = size (ur);
	if (nc != 1)
		ur = ur.';
	endif
	U = L \ ur;
	U = U.';	% let's keep all our vectors as row-vectors

	% we still have to add the homogenous solution
	% first we add the constant so U(1) = 0
	n = U(1);
	U = U .- n;

	% now the linear part
	k = (1 - U(N))/rmax;

	U = U + k*r;

	return;
endfunction

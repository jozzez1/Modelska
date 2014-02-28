function [E, u, U] = iterate (E0, D, rmin, rmax, N, prec)
	% initial guess
	[E1, u1, U1] = ionStart (rmin, rmax, N);

	% first step
	[E2, u2, U2] = step (u1, U1, E1, D, rmin, rmax, N);

	k = 0;
	while (sumsq(u1 - u2) + sumsq(U1 - U2) > prec)
		E1 = E2;
		u1 = u2;
		U1 = U2;
		[E2, u2, U2] = step (u1, U1, E0, D, rmin, rmax, N);
		k++;

		if k > 1e3
			printf ("Failed to converge in 1000 iterations.\n");
			break;
		endif
		k,
		sumsq(u1 - u2) + sumsq(U1 - U2)
	endwhile

	E = E2;
	u = u2;
	U = U2;

	return;
endfunction

function u = correct (ell, rmin, rmax, N)
	x = linspace (rmin, rmax, N);
	if ell == 0
		u = x .* exp (-x);
	elseif ell == 1
		u = (x .^ 2) .* exp (-0.5 .* x);
	elseif ell == 3
		u = x .* (1 - 0.5 * x) .* exp (-0.5 .* x);
	endif

	norma = sqrt(sumsq(u));
	u = (u ./ norma) * sqrt(N / (rmax - rmin));

	return;
endfunction

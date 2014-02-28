function E = energy (R, phi, rmin, rmax, N)
	r = linspace (rmin, rmax, N);
	h = (rmax - rmin)/N;
	Z = 2;
	E0 = -13.6;

	derivR(1) = 0;
	for i = 2:N
		derivR (i) = (R(i) - R(i-1))/h;
	end

	f = derivR.^2 - 2*Z*((R.^2) ./ r) - (R.^2 .* (phi ./ r));
	f(1) = 0;

	E = 2 * E0 * sum(f) * h;
	return;
endfunction

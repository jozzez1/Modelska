function z = solve1 (N, md)
	pkg load signal
	A = dens (N, md);

	%% for now let's compute 2D Four. trans
	Fhat = dst2 (A);
	
	for m = 1:N
		for n = 1:N
			zhat (m,n) = (1/N)^2 * Fhat (m,n);
			zhat (m,n) /= 2 * (cos(pi*m/N) + cos(pi*n/N) - 2);
		end
	end

	z = idst2 (zhat);

	pkg unload signal
	return;
endfunction

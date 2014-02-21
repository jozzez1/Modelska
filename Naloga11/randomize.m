function M = randomize (A)
	N = size(A)(1);
	D = sum(sum(A))/N - 1;
	Imax = N*D;

	for i = 1:Imax
		M = randcon (A);
		A = M;
	end

	return;
endfunction

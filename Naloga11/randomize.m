function M = rand_con (N, D)
	M = begin (N, D);
	Imax = (N*D/2)^2;

	for i = 1:Imax
		M = randomize (M);
	end
endfunction

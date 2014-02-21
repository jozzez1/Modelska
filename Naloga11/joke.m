function [t, r] = joke (M)
	N = size (M)(1);
	tmax = radiusz (M) * N;
	tmax = 10;

	x = zeros (1,N);
	x (1) = 1;	% joke starts from the "patient" zero
	y = x;

	t (1) = 0;
	r (1) = sum(x)/N;
	p = 0.2;	% probability that the joke will spread

	for k = 2:tmax
		for j = 1:N
			% someone knows the "joke" and he will tell it to his friends!
			if x(j) == 1
				for i = 1:N
					if (M(i,j) == 1) && (p >= rand())
						y(i) = 1;
					endif
				end
			endif
		end
		x = y;
		t(k) = k-1;
		r(k) = sum(x)/N;
	end

	return;
endfunction

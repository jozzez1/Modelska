function M = begin (N, D)
	assert (mod(D,2) == 0);

	M = eye (N); % we need this, or it will be wrong
%	M = zeros(N);
	% now we add the diagonals
	for i = 1:D/2
		M += diag(ones(1, N-i), +i);
		M += diag(ones(1, N-i), -i);
	end

	% now we have to add "islands" on edges, due to the periodicity
	for i = 1:D/2
		M += diag(ones(1, i), N-i);
		M += diag(ones(1, i), i-N);
	end

	return;
endfunction

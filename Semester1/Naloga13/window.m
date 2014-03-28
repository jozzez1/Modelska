function H = window (N, M)
	H = 1 - hamming (N);

	for i = 1:N
		if (i > N/2-M) && (i < N/2+M)
			H(i) = 0;
		endif
	end

	return;
endfunction

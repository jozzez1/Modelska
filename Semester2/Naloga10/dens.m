function A = dens (N, md)
	%% this function returns the density distribution,
	%% depending on which "mode" we choose

	A = [];
	if md != 4
		for i = 1:N
			for j = 1:N
				A(i,j) = 1;
				if (i < 0.5*N) && (i > 0.25*N) && (j > 0.5*N) && (j < 0.75*N)
					A (i,j) += md;
				endif
			end
		end
	elseif md == 4		%% cylinder boundaries
		A = zeros (N,N);
		A (1,:) = -1/(N^2);
		A (N,:) = -1/(N^2);
	endif

	return;
endfunction

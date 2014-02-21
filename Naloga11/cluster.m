function C = cluster (M)
	N = size(M)(1);
	C = 0;

	% so we start
	for i = 1:N		% column
		for j = 1:N
			if j == i
				continue;
			endif
			for k = 1:N
				if (k == i) || (k == j)
					continue;
				endif
				
				if (M (j, i) == 1) && (M (k, i) == 1) && (M (k, j) == 1)
					C++;
				endif
			end
		end
	end

	% we also counted all the permutations,
	% so we have to divide with 3! = 6
	C /= 6;

	return;
endfunction

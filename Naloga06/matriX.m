function [X, j] = matriX (T, k)
	% we take the good selected parts of T and write it to X
	n = dec2bin (k,10);

	j = 1;
	X = [];
	for i = 1:10
		if bin2dec(n(i)) == 1
			X(:,j) = T(:,i);
			j++;
		endif
	end

	j = j-1;

	return;
endfunction

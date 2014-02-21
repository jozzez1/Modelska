function [chi2, k, S] = allModes (T, R, y)
	for i = 1:1023
		[b, chi2(i), P, S{i}, k(i), n] = regression (T, R, y, i);
	end
	return;
endfunction

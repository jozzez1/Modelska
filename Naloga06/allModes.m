function [chi2, k] = allModes (T, R, y)
	for i = 1:1023
		[b, chi2(i), P, X, k(i), n] = regression (T, R, y, i);
		i
	end
	return;
endfunction

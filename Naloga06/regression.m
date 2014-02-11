function [b, chi2, P, S, j, n] = regression (T, R, y, k)
	[X, j] = matriX (T, k);
	[U, S, V] = svd (X, 1);
	n = dec2bin (k, 10);

	P = inv(transpose(X) * inv(R) * X);
	b = P * transpose(X) * inv(R) * y;

	chi2 = (y - X*b).' * inv(R) * (y - X*b);

	chi2 /= (10 - j - 1);

	return;
endfunction

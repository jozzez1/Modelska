function plotter2 (conf, plotflag, model)
	[T, R, y] = matrixTwr (model);

	[b, chi2, P, S, j, n] = regression (T, R, y, conf)

	beta = zeros (1,10);
	k = 1;
	for i = 1:10
		beta(i) = bin2dec (n(i));
		if beta(i) == 1
			beta(i) *= b(k);
			k++;
		endif
	end

	x = linspace (-1.5, 1.5, 100);
	y = beta(1) .+ x*beta(2) + (x .^ 2)*beta(3) + (x .^ 3)*beta(4) + (x .^ 4)*beta(5) + (x .^ 5)*beta(6) + (x .^ 6)*beta(7) + (x .^ 7)*beta(8) + (x .^ 8)*beta(9) + (x .^ 9)*beta(10);

	plot (x, y);

	if plotflag == 1
		M(:,1) = x.';
		M(:,2) = y.';

		tit = sprintf("%d-%d.txt", conf, model);
		save(tit, "-ascii", "M");
	endif
endfunction

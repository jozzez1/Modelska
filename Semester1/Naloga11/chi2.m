function val = chi2 (args)
	% args is a vector of data
	M = load ('-ascii', 'zgodovina1.dat');
	T = M(:,1);
	N = M(:,2);

	y = model (T, args);
	v = y - N;
	%K, tau, T0

	val = v.' * v/(length(T) - 4);
	
	return;
endfunction

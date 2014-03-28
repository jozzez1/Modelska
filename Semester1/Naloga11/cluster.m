function C = cluster (M)
	% this is actually just a wrapper for a program I wrote in C
	% sources are in cluster.c + Makefile
	% DO NOT use -O3 to compile!

	% first write it to a file
	mprint (M);

	% then the C program will process it
	N = size (M)(1);
	command = sprintf ("./ccluster %d", N);
	[out, result] = system (command);
	system ("rm *.txt");

	% and now we convert it back to number
	C = str2num (result);

	return;
endfunction

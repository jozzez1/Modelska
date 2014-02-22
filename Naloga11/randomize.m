function M = randomize (A, save_progress)
	% wrapper for crandomize (source is randomize.c) function, since it's way faster
	N = size(A)(1);
	D = sum(sum(A))/N - 1;
	
	mprint (A);
	command = sprintf("./crandomize %d %d", N, save_progress);
	system (command);
	mfile = sprintf('matrix-%d.txt', N);
	M = load ('-ascii', mfile);

	command = sprintf ("rm %s", mfile);
	system (command);

	return;
endfunction

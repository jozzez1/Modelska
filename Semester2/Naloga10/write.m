function write (N, md)
	%% to save
	z = solve1 (N, md);

	%% we add the boundaries
	if (md != 4)
		Z = horzcat (zeros (N,1), z, zeros (N,1));
		Z = horzcat (zeros (N+2,1), Z', zeros (N+2, 1))';
	elseif (md == 4)
		Z = horzcat (zeros(N,1), z, zeros(N,1));
		Z = horzcat (ones(N+2,1), Z', ones(N+2,1))';
	endif
	h = 1/N;
	X = horzcat(0, linspace (h/2, 1-h/2, N), 1);
	Y = X;

	[XX, YY] = meshgrid (X, Y);

	surf (XX, YY, Z);
	view (2);

	XX = reshape (XX, (N+2)*(N+2), 1);
	YY = reshape (YY, (N+2)*(N+2), 1);
	ZZ = reshape (Z, (N+2)*(N+2), 1);

	A = [XX, YY, ZZ];
	data = ["solution-", 'N', "-", 'md', ".txt"]

	save ('-ascii', data, 'A');
endfunction

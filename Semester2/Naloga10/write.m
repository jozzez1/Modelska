function write (N, md, hw)
	%% to save
    base = "solution";
    switch (hw)
        case 1
        	z = solve1 (N, md);
        case 2
            z = solve2 (N, md);
        case 3
            z = solve3 (N, md);
            base = "cilinder";
    endswitch

	%% we add the boundaries
	if (md != 4)
		Z = horzcat (zeros (N,1), z, zeros (N,1));
		Z = horzcat (zeros (N+2,1), Z', zeros (N+2, 1))';
	elseif (md == 4)
		Z = horzcat (zeros(N,1), z, zeros(N,1));
		Z = horzcat (ones(N+2,1), Z', ones(N+2,1))';
	endif
	h = 1/N;
	Y = horzcat(0, linspace (h/2, 1-h/2, N), 1);
    if hw != 3
    	X = Y;
    elseif hw == 3
        X = horzcat(-0.5, linspace (-0.5+h/2, 0.5-h/2, N), 0.5);
    endif

	[XX, YY] = meshgrid (X, Y);
	XX = reshape (XX, (N+2)*(N+2), 1);
	YY = reshape (YY, (N+2)*(N+2), 1);
	ZZ = reshape (Z, (N+2)*(N+2), 1);

	A = [XX, YY, ZZ];
	data = ["solution_", int2str(N), "_", num2str(md), ".txt"]
	save ('-ascii', data, 'A');

    system (["./prg ", data, " ",int2str(hw)]);
endfunction

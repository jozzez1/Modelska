function M = randomize (A)
	N = size(A)(1);
	% we randomly select one of the nodes (column)
	j = randi(N);

	% now we have to select one of it's connections (row)
	i = randi(N);
	while (A(i, j) == 0) || (i == j)
		i = randi(N);
	endwhile

	% now we select place to move it to
	k = randi(N);
	while (A(k,j) == 1) || (k == j)
		k = randi(N);
	endwhile

	% and now we make the move
	% 1.) we cut the former connection j <--> i
	A (i,j) = 0;
	A (j,i) = 0;
	% 2.) we assign a new connection j <--> k
	A (k,j) = 1;
	A (j,k) = 1;

	M = A;
	return;
endfunction

function M = randcon (A)
	N = size(A)(1);

	j = randi (N);
	while (sum(A(:,j)) == N) || (sum(A(:,j)) == 1)
		j = randi (N);
	end

	% to prevent dead loops we will first
	% check which connections exist
	x = [];
	y = [];
	kx = 1;
	ky = 1;
	for i = 1:N
		if (A(i,j) == 1) && (i != j)
			x(kx) = i;
			kx++;
		elseif (A(i,j) == 0)
			y(ky) = i;
			ky++;
		endif
	end

	% now we pick one of the connections to move
	Nx = kx-1;
	Ny = ky-1;

	i = x(randi (Nx));
	k = y(randi (Ny));

	% and now we move the stuff
	A (i,j) = 0;
	A (j,i) = 0;

	A (k,j) = 1;
	A (j,k) = 1;

	M = A;
	return;
endfunction

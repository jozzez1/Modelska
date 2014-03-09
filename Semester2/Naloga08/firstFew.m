function E = firstFew (N)
	[A, B] = matrixAB (1, N);
	e0 = eigs (A, B, 3, "sm");
	
	[A, B] = matrixAB (2, N);
	e1 = eigs (A, B, 3, "sm");
	
	[A, B] = matrixAB (3, N);
	e2 = eigs (A, B, 3, "sm");

	[A, B] = matrixAB (4, N);
	e3 = eigs (A, B, 3, "sm");

	[A, B] = matrixAB (5, N);
	e4 = eigs (A, B, 3, "sm");

	E = [e0; e1; e2; e3; e4];
	[E, i] = sort (E);

	E = E(1:13,:);

	return;
endfunction

function predict (filename, order)
	warning ("off", "Octave:shadowed-function");
	pkg load tsa;
	warning ("on", "Octave:shadowed-function");

	y = load ('-ascii', filename);
	N = length(y);
	
	% we take only half of the data
	t = y(N/2+1:N,:);
	y = y(1:N/2,:);
	N = length(y);

	% get the coefficients
	[A, RC, PE] = lattice (y.', order);
	sigma2 = max(PE);

	for k = 1:N
		y(k+N) = 0;
		for j = 1:order
			y(k+N) += A(j) * y(N + k - j);
		end
	end

	y = y(N+1:2*N,:);

	x = 0:N-1;
	err = sqrt(sigma2)*ones(N,1);

	Result = [x.', t, y, err];
	output = [filename, "-predict.dat"];

	save ('-ascii', output, 'Result');

	return;
endfunction


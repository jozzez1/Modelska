function [y, Py] = filtered (filename, err)
	y = load ('-ascii', filename);

	fy = fft (y);

	% Wiener filter
	W = 1 ./ (1 + err ./ (abs(fy - sqrt(err)).^2));

	% and now the spectrum
	Py = W .* (abs(fy).^2);

	return;
endfunction

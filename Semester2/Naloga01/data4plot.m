function data4plot (filename)
	% let's get the approximations for a few selections
	[P10, z10, A10, sigma210] = spect (filename, 10);
	[P20, z20, A20, sigma220] = spect (filename, 20);
	[P30, z30, A30, sigma230] = spect (filename, 30);
	[P40, z40, A40, sigma240] = spect (filename, 40);

	% and now we'll use the Wiener filter with the best error estimate
	[y, Py] = filtered (filename, sigma240);
	w = linspace(0, 2*pi*(length(Py)-1)/length(Py), length(Py));

	Result = [w.', P10.', P20.', P30.', P40.', Py];
	
	output = [filename, "-spectra.dat"];
	save ('-ascii', output, 'Result');
endfunction

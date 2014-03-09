function triplotter (a)
	M = load ('-ascii', ['mode-', a, '.dat']);
	t = load ('-ascii', 'triangles.dat');
	X = M (:,1);
	Y = M (:,2);
	Z = M (:,3);
	trisurf (t, X, Y, Z);
	view(2);

	orient ("landscape");
	axis ("off");
	print (['mode-', a, '.png'], '-S1200,600');
endfunction

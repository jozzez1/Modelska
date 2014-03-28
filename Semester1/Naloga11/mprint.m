function mprint (M)
	N = size (M)(1);
	dat = sprintf ("matrix-%d.txt", N);

	fout = fopen (dat, "w");
	for i = 1:N
		fprintf (fout, "%d ", M(i,:));
		fprintf (fout, "\n");
	end
	fclose (fout);
endfunction

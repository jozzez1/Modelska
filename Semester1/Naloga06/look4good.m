function [n, m] = look4good (chi2, k, S)
	j = 1;

	fout = fopen ("good2check-eight.txt", "w");
	for i = 1:1023
		if (chi2(i) < 1.4) && (chi2(i) > 0.6) && (k(i) < 7) && (diag(S{i})(k(i)) > 1)
			n(j,1) = i;
			n(j,2) = chi2(i);
			n(j,3) = k(i);
			n(j,4) = diag(S{i})(k(i));

			m = dec2bin(i,10);

			fprintf (fout, "%d  %e  %e  %e  %s\n",
				n(j,1), n(j,2), n(j,3), n(j,4), m);

			j++;
		endif
	end

	%save ("good2check-theta.txt", '-ascii', "n", "m");
	fclose (fout);

	return;
endfunction

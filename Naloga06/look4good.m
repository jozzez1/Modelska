function n = look4good (chi2, k, S)
	j = 1;
	for i = 1:1023
		if (chi2(i) < 8) && (chi2(i) > 0.5) && (k(i) < 7) && (diag(S{i})(k(i)) > 1)
			n(j,1) = i;
			n(j,2) = chi2(i);
			n(j,3) = k(i);
			n(j,4) = diag(S{i})(k(i));
			j++;
		endif
	end

	save -ascii good2check.txt n

	return;
endfunction

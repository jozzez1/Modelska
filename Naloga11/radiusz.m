function r = radiusz (M)
	r = 1;
	m = min (min(M));
	while m == 0
		r++;
		m = min (min (M^r));
	endwhile

	return;
endfunction

function F = idst2 (Fhat)
	F = idst(idst(Fhat)')';
	return;
endfunction

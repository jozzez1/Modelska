function [T, R, y] = matrixTwr (k)
	M = [];
	if k == 1
		M = load ("wrboost_eight.dat", '-ascii');
	elseif k == 2
		M = load ("wrboost_omega.dat", '-ascii');
	elseif k == 3
		M = load ("wrboost_theta.dat", '-ascii');
	endif

	for i = 1:10
		T(:,i) = M(:,1) .^ (i-1);
	end

	[nr, nc] = size (T);
	R = 2 * eye(nr);

	y = M(:,2);

	return;
endfunction

function [A, B] = matrixAB (m, N)
	for i = 1:N
		for j = 1:N
			A(i,j)  = (i*j)/(2*m + i + j);
			A(i,j) -= (2*i*j + i + j)/(2*m + i + j + 1);
			A(i,j) += ((i + 1)*(j + 1))/(2*m + i + j + 2);

			B(i,j)  = 1/(2*m + i + j + 2);
			B(i,j) -= 2/(2*m + i + j + 3);
			B(i,j) += 1/(2*m + i + j + 4);
		end
	end

	return;
endfunction

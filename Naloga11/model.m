function val = model (x, args)
	K = args(1);
	tau = args(2);
	T0 = args(3);

	val = (K^2) * atan (tan((x .- T0)./(K * tau) .+ pi/2)./K);

	N = length (val);
	for i = 1:N
		if (x(i) > T0)
			val(i) += K^2 * pi;
		endif
	end

	return;
endfunction


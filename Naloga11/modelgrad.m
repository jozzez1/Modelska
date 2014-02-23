function dNdp = modelgrad (x,f,p,dp,func)
	K = p(1);
	tau = p(2);
	T0 = p(3);

	prod1 = K^2 ./ (K^2 + (tan ((x - T0)/(K*tau) - pi/2)).^2);
	cosstuff = 1 ./ (cos((x - T0)./(K*tau) - pi/2).^2);
	
	dNdK = - prod1 .* (((x - T0)./(K*tau)).*cosstuff .+ tan((x .- T0)./(K*tau) - pi/2));
	dNdtau = - prod1 .* ((x - T0)./(tau^2)) .* cosstuff;
	dNdT0 = - (prod1 ./ tau) .* cosstuff;

	for i = 1:length (dNdK)
		if (x(i) > T0)
			dNdK (i) += K^2 * pi;
		endif
	end

	if (isrow(dNdK))
		dNdK = dNdK';
	endif

	if (isrow(dNdtau))
		dNdtau = dNdtau';
	endif

	if (isrow(dNdT0))
		dNdT0 = dNdT0';
	endif

	dNdp = [dNdK, dNdtau, dNdT0];
	return;
endfunction

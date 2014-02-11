function plotter (conf)

[T, R, y] = matrixT ();
[b, chi2, P, S, j, n] = regression (T, R, y, conf)

beta = zeros (1,10);
k = 1;

for i = 1:10
	beta(i) = bin2dec (n(i));
	if beta(i) == 1
		beta(i) *= b(k);
		k++;
	endif
end

T = linspace (90, 362, 41);
P = linspace (272, 602, 41);

[TT, PP] = meshgrid (T, P);

z = beta(1) .+ TT*beta(2) + PP*beta(3) + (TT.^2)*beta(4) + (TT .* PP)*beta(5) + (PP .^ 2)*beta(6) + (TT .^ 3)*beta(7) + ((TT .^ 2) .* PP)*beta(8) + (TT .* (PP .^ 2))*beta(9) + (PP .^ 3)*beta(10);

meshz (T, P, z);
view (130, 30);

xlabel ("Temperatura -- $T$");
ylabel ("Elektricna moc -- $P$");
zlabel ("Toplotna prevodnost -- $\\lambda$");
tit = sprintf ("Rezultat za konfiguracijo %d, tj. %s", conf, n);
title (tit);

% we add the points on which we fit
hold on;

plot3 (9.00000000e+01, 2.76000000e+02, 4.23450000e+01, "+");
plot3 (1.00000000e+02, 5.45000000e+02, 4.16000000e+01, "+");
plot3 (1.49000000e+02, 2.75000000e+02, 3.95375000e+01, "+");
plot3 (1.61000000e+02, 6.02000000e+02, 3.77875000e+01, "+");
plot3 (2.06000000e+02, 2.74000000e+02, 3.73525000e+01, "+");
plot3 (2.27000000e+02, 5.38000000e+02, 3.64975000e+01, "+");
plot3 (2.47000000e+02, 2.74000000e+02, 3.63600000e+01, "+");
plot3 (2.70000000e+02, 5.50000000e+02, 3.57850000e+01, "+");
plot3 (3.52000000e+02, 2.72000000e+02, 3.39150000e+01, "+");
plot3 (3.62000000e+02, 5.22000000e+02, 3.45300000e+01, "+");

axis ([50 370 200 610]);

hold off;

endfunction


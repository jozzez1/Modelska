reset

set xlabel 'bredzim. radij -- $x$'
set ylabel 'radialna funkcija -- $R(x)$'


set term epslatex color solid size 12cm,18cm

set output "radials.tex"
set multiplot layout 2,1
set xrange [0:20]
set title 'Rezultati radialnih funkcij, $R(x)$'
plot "radial0.dat" u 1:2 w l title '$n = 1,\ \ell = 0$', \
"radial20.dat" u 1:2 w l title '$n = 2,\ \ell = 0$', \
"radial21.dat" u 1:2 w l title '$n = 2,\ \ell = 1$'

set ylabel 'potencial -- $\varphi(x)$'
set xrange [0:10]
set title 'Rezultati potencialov, $\varphi(x)$'
plot "radial0.dat" u 1:4 w l title '$n = 1,\ \ell = 0$', \
"radial20.dat" u 1:4 w l title '$n = 2,\ \ell = 0$', \
"radial21.dat" u 1:4 w l title '$n = 2,\ \ell = 1$'
unset multiplot
unset output

set output "errors.tex"
set multiplot layout 2,1
set autoscale
set title 'Logaritem absolutne napake radialnih funkcij'
set ylabel '$\log_{10}|R(x) - \hat{R}(x)|$'
plot "radial0.dat" u 1:(log10(abs($2-$3))) w l title '$n = 1,\ \ell = 0$', \
"radial20.dat" u 1:(log10(abs($2-$3))) w l title '$n = 2,\ \ell = 0$', \
"radial21.dat" u 1:(log10(abs($2-$3))) w l title '$n = 2,\ \ell = 1$'

set ylabel '$\log_{10}|\varphi(x) - \hat{\varphi}(x)|$
set title 'Logaritem absolutne napake potencialov'
plot "radial0.dat" u 1:(log10(abs($4-$5))) w l title '$n = 1,\ \ell = 0$', \
"radial20.dat" u 1:(log10(abs($4-$5))) w l title '$n = 2,\ \ell = 0$', \
"radial21.dat" u 1:(log10(abs($4-$5))) w l title '$n = 2,\ \ell = 1$'
unset multiplot
unset output

set output "radials.tex"
set multiplot layout 2,1
set xrange [0:20]
set title 'Rezultati radialnih funkcij, $R(x)$'
plot "radial0.dat" u 1:2 w l title '$n = 1,\ \ell = 0$', \
"radial20.dat" u 1:2 w l title '$n = 2,\ \ell = 0$', \
"radial21.dat" u 1:2 w l title '$n = 2,\ \ell = 1$'

set ylabel 'potencial -- $\varphi(x)$'
set xrange [0:10]
set title 'Rezultati potencialov, $\varphi(x)$'
plot "radial0.dat" u 1:4 w l title '$n = 1,\ \ell = 0$', \
"radial20.dat" u 1:4 w l title '$n = 2,\ \ell = 0$', \
"radial21.dat" u 1:4 w l title '$n = 2,\ \ell = 1$'
unset multiplot
unset output

set term epslatex color solid size 12cm,9cm
set output "helij.tex"
set autoscale
set xlabel 'brezdim. radij $x$'
set ylabel 'Amplituda'
set title 'Helijev atom'

plot "helium.dat" u 1:2 w l title '$R(x)$', \
"helium.dat" u 1:(-$3) w l lt 3 title '$-\varphi(x)$'
unset output


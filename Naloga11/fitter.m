M1 = load ('-ascii', 'mestaSLO.txt');
M2 = load ('-ascii', 'mestaZDA.txt');
M3 = load ('-ascii', 'mestaIND.txt');
M4 = load ('-ascii', 'mestaCHN.txt');

N1 = size(M1)(1);
N2 = size(M2)(1);
N3 = size(M3)(1);
N4 = size(M4)(1);

X1(:,1) = ones (N1,1);
X2(:,1) = ones (N2,1);
X3(:,1) = ones (N3,1);
X4(:,1) = ones (N4,1);

X1(:,2) = log(M1(:,1));
X2(:,2) = log(M2(:,1));
X3(:,2) = log(M3(:,1));
X4(:,2) = log(M4(:,1));

y1 = log(M1(:,2));
y2 = log(M2(:,2));
y3 = log(M3(:,2));
y4 = log(M4(:,2));

P1 = inv (X1.' * X1);
P2 = inv (X2.' * X2);
P3 = inv (X3.' * X3);
P4 = inv (X4.' * X4);

% calculated parameters
b1 = P1 * X1.' * y1;
b2 = P2 * X2.' * y2;
b3 = P3 * X3.' * y3;
b4 = P4 * X4.' * y4;

% relative errors
err1 = abs(diag(sqrt(P1)) ./ b1);
err2 = abs(diag(sqrt(P2)) ./ b2);
err3 = abs(diag(sqrt(P3)) ./ b3);
err4 = abs(diag(sqrt(P4)) ./ b4);

%x = 1:280;
%f1 = exp(b1(1)) .* x.^(b1(2));
%f2 = exp(b2(1)) .* x.^(b2(2));
%f3 = exp(b3(1)) .* x.^(b3(2));
%f4 = exp(b4(1)) .* x.^(b4(2));

%A = [x', f1', f2', f3', f4'];

%save -ascii cities.txt A
f1 = X1 * b1;
f2 = X2 * b2;
f3 = X3 * b3;
f4 = X4 * b4;

chi1 = (f1 - y1).' * (f1 - y1);
chi2 = (f2 - y2).' * (f2 - y2);
chi3 = (f3 - y3).' * (f3 - y3);
chi4 = (f4 - y4).' * (f4 - y4);

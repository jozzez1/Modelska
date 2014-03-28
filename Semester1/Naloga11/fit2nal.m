% our test function
% data
M = load ('-ascii', 'zgodovina.dat');
N = M(:,2);
T = M(:,1);

% model function is model.m
% initial values
init = [5000; 100; 2020];

F = @model;
%dFdp = @modelgrad; -- doesn't work well
dFdp = 'dfdp';	% we will use numerical derivatives

wt = ones (length(N),1);

global verbose;
verbose = 1;

niter = 200;
stol = 0.0001;
dp = [1e-5;1e-5;1e-5];

[f, p, cvg, iter, corp, covp, covr, stdresid, Z, r2] = leasqr (T, N, init, F, stol, niter, wt, dp, dFdp)

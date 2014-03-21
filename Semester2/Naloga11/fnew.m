function fT = fnew (FT, f1, f0, ht)
    N = length(f0);
    h = 1.0 / N;

    for i = 2:N-1
        fT(i) = 0.5 * (f1(i+1) - f1(i-1)) * (FT(i+1) - FT(i-1)) + FT(i)*(f1(i+1) + f1(i-1) - 2*f1(i));
        fT(i) = fT(i) * ((ht/h)^2) + 2*f1(i) - f0(i);
    end

    fT(1) = 2*cos(fT(2))/(N * FT(2)) - fT(3);
    fT(N) = 2*fT(N-1) - fT(N-1);

    %% we have to fix it to become a column-vector
    fT = fT';

    return;
endfunction

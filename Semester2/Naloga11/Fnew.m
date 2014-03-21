function FT = Fnew (f1, f0, ht)
    N = length(f0);
    h = 1.0 / N;

    A = diag(ones(N-3,1),1) + diag(ones(N-3),-1);
    for i = 2:N-1
        v(i-1) = 2 + 0.25*(f1(i+1) - f1(i-1))^2;
        df(i-1) = ((h/ht)^2) * (f1(i-1) - f0(i-1));
    end
    A -= diag(v);
    FT = A \ df;
    
    %% we fix the borders
    F0 = h*sin(f1(1)) - FT(1);
    FN = 0;

    FT = [F0; FT; FN];

    return;
endfunction

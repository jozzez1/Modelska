function FT = Fnew (f1, f0, ht)
    N = length(f0);
    h = 1.0 / N;

    A = diag(ones(N-3,1),1) + diag(ones(N-3, 1),-1);
    v(1) = 1 + 0.25*(f1(3) - f1(1))^2;
    df(1) = ((h/ht)^2) * (f0(2) - f1(2))^2 - h*sin(f1(1));
    for i = 3:N-1
        v(i-1) = 2 + 0.25*(f1(i+1) - f1(i-1))^2;
        df(i-1) = ((h/ht)^2) * (f0(i) - f1(i))^2;
    end
    df = df';
    A -= diag(v);
    FT = A \ df;
    
    %% we fix the borders
    F0 = h*sin(f1(1)) - FT(1);
    FN = 0;

    FT = [F0; FT; FN];

    return;
endfunction

function z = solve3 (N, md)
    pkg load signal

    F = dens (N, 4);
    Fphat = dst (F).';

    %% ok, time for phun
    h = 1.0/N;
    r = linspace (-0.5+h/2,0.5-h/2, N);
    a = diag (ones(1,N-1) - 0.5 * h ./ r(1,2:N), -1) + diag (ones(1,N-1) + 0.5 * h ./ r(1,1:N-1), +1);

    for k = 1:N
        A = a + diag (ones(N,1) * (2 * cos(pi*k/N) - 4));
        zhat (:,k) = A \ (Fphat (:,k) ./ (N^2));
    end

    z = idst (zhat.');
    return;
endfunction

function z = solve2 (N, md)
    pkg load signal

    F = dens (N, md);
    Fphat = dst (F).';

    %% now we take the colums and save a system of equations
    %% A * U = F, where U, F and A are N \times N matrices
    %% so let's construct matrix A
    a = diag (ones(N-1, 1), +1) + diag(ones(N-1, 1), -1);
    
    for k = 1:N
        A = a + diag (ones(N,1) * (2*cos(pi*k/N) - 4));
        zhat (:,k) = A \ (Fphat (:,k) ./ (N^2));
    end
    
    z = idst (zhat.');

    pkg unload signal
    return;
endfunction

function [FF, ff] = solve (N, ht, fi0, T)
    f0 = start (fi0, N);
    f1 = f0;

    ff = [f0, f1];
    FF = zeros (N,1);
    FT = Fnew (f1, f0, ht);

    FF = [FF, FT];

    for i = 1:T
        F = FT;
        [fT, FT] = time_step (F, f1, f0, ht);

        f0 = f1;
        f1 = fT;

        ff = [ff, fT];
        FF = [FF, FT];
    end

    return;
endfunction

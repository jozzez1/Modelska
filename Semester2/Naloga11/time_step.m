function [fT, FT] = time_step (F, f1, f0, ht)
    fT = fnew (F, f1, f0, ht);
    FT = Fnew (fT, f1, ht);

    return;
endfunction

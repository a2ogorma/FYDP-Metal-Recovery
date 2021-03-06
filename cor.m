function func = cor(Ecorr, Erev, iLa, iLc, S_PCB, on_PCB_cat, on_PCB_an, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alphas z
    i_BV_corr = i_BV(Ecorr-Erev, i0, iLa, iLc, alphas, z, temp);
    i_corr_an = on_PCB_an.*subplus(i_BV_corr);
    i_corr_cat = on_PCB_cat.*(-subplus(-i_BV_corr));
    i_corr = i_corr_an + i_corr_cat;
    %arrange surface areas in proper order
    S_corr = [S_PCB(2:3) sum(S_PCB(2:7)) S_PCB(4:5) S_PCB(5:6) S_PCB(6:7) sum(S_PCB(2:7)) sum(S_PCB(2:7))];
    if abs(sum(S_corr.*i_corr)) == Inf
        disp(S_corr)
        disp(i_corr)
        error("Corrosion currents infinite");
    end
    func = sum(S_corr.*i_corr);
end
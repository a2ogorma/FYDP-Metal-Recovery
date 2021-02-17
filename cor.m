function func = cor(Ecorr, Erev, iL, S_PCB, on_PCB_cat, on_PCB_an, psbl_lch, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alphas z
    i_BV_corr = i_BV(Ecorr-Erev, i0, iL, alphas, z, temp);
    i_corr_an = on_PCB_an.*subplus(i_BV_corr);
    i_corr_cat = psbl_lch.*on_PCB_cat.*(-subplus(-i_BV_corr));
    i_corr = i_corr_an + i_corr_cat;
    %arrange surface areas in proper order
    S_corr = [S_PCB(2:5) sum(S_PCB(2:9)) S_PCB(6:9) sum(S_PCB(2:9)) sum(S_PCB(2:9))];
    if abs(sum(S_corr.*i_corr)) == Inf
        disp(S_corr)
        disp(i_corr)
        error("Corrosion currents infinite");
    end
    func = sum(S_corr.*i_corr);
end
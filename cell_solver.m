function Y = cell_solver(I, E_an, E_cat, V_app, r_sol, r_hardware, Erev, iL, onCathode, onAnode, psbl_cell, S_an, S_cat, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    %onCathode: Array of reaction locations based on Erev. 1
    %corresponds to reaction occurring on cathode, 0 corresponds to
    %reaction occuring on anode
    global i0 alphas z
    eta_cat = E_cat - Erev;
    eta_an = E_an - Erev;
    i_cat_temp = onCathode.*(i_BV(eta_cat, i0, iL, alphas, z, temp));
    i_cat_cat = psbl_cell(1,:).*-subplus(-i_cat_temp);
    i_cat_an = psbl_cell(2,:).*subplus(i_cat_temp);
    i_cat = i_cat_cat + i_cat_an;
    I_cat = i_cat*S_cat;
    if sum(I_cat) == -Inf
        error("Cathodic current infinite");
    end
    i_an = onAnode.*i_BV(eta_an, i0, iL, alphas, z, temp);
    I_an = i_an*S_an;
    if sum(I_an) == Inf
        error("Anodic current infinite");
    end
    Y(1) = E_an - E_cat + I*(r_sol+r_hardware) - V_app;
    Y(2) = I - sum(I_an);
    Y(3) = sum(I_an) + sum(I_cat);
end
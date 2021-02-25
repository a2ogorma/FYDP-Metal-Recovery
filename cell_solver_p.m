function Y = cell_solver_p(I, E_an, E_cat, V_app, r_sol, r_hardware, Erev_cat, Erev_an, iLa_cat, iLc_cat, iLa_an, iLc_an, onCathode, onAnode, S_an, S_cat_p, temp)
    %potentiostat
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    %onCathode: Array of reaction locations based on Erev. 1
    %corresponds to reaction occurring on cathode, 0 corresponds to
    %reaction occuring on anode
    global i0 alphas z
    S_cat = sum(S_cat_p);
    eta_cat = E_cat - Erev_cat;
    eta_an = E_an - Erev_an;
    i_cat = onCathode.*(i_BV(eta_cat, i0, iLa_cat, iLc_cat, alphas, z, temp));
    i_cat_cat = -subplus(-i_cat);
    i_cat_an = subplus(i_cat);
    I_cat_cat = i_cat_cat*S_cat;
    I_cat_an = i_cat_an.*[S_cat_p(1:2) S_cat S_cat_p(3:4) S_cat_p(4:5) S_cat_p(5:6) 0 S_cat];
    I_cat = I_cat_cat+I_cat_an;
    if sum(I_cat) == -Inf
        error("Cathodic current infinite");
    end
    i_an = onAnode.*i_BV(eta_an, i0, iLa_an, iLc_an, alphas, z, temp);
    I_an = i_an*S_an;
    if sum(I_an) == Inf
        error("Anodic current infinite");
    end
    Y(1) = (E_an - E_cat + I*(r_sol+r_hardware) - V_app)*1E6;
    Y(2) = (I - sum(I_an))*1E6;
    Y(3) = (sum(I_an) + sum(I_cat))*1E6;
end
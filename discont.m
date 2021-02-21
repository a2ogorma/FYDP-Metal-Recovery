function [flag, isterminal, direction] = discont(t, Cm, temp, pres, vol_cell, vol_lch, Q, S_an, S_cat, mode, VI_app, n_particles, l, A_cell)
    direction = [];
    isterminal = zeros(1,13);
    %isterminal(12) = 1;
    %isterminal(13) = 1;
    flag = ones(1,13);
    global F z i0 km alphas lamda rho aH2 aO2
    if mode == 1 %potentiostat
        V_app = VI_app;
    elseif mode == 2 %galvanostat
        I_app = VI_app;
    end
    %calculate total mass and wt fractions based on partial masses at
    %timestep
    m_PCB = (Cm(25:33).');
    m_PCB_total = sum(m_PCB);
    %Convert weight to volume and volume fraction
    V_PCB = m_PCB./rho; %Element volumes in m^3
    vfrac_PCB = V_PCB/sum(V_PCB);
    V_PCB_total = sum(V_PCB);
    
    %calculate particle radius in m
    r_particles = (V_PCB_total/(n_particles*4*pi)*3)^(1/3);
    
    %Surface area calculation for corrosion
    packing_density = 0.6;
    SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
    S_PCB_total = V_PCB_total*packing_density*SSA; %m2 Total surface area of spheres in lch.
    S_PCB = vfrac_PCB*S_PCB_total; %array with exposed surface area of each component
    S_PCB = subplus(S_PCB*100^2); %convert m^2 to cm^2
    
    [Erev_cell, Erev_lch, psbl_cell, psbl_lch] = nernstPotential(Cm,temp,aH2,aO2);
    
    %resistance calculation for IR drop in cell
    kappa = 1000*(Cm(1)*lamda(1)+Cm(2)*lamda(2)+Cm(3)*lamda(3)+...
        Cm(4)*lamda(4)+Cm(5)*lamda(5)+Cm(6)*lamda(6)+Cm(7)*lamda(7)+...
        Cm(8)*lamda(8)+Cm(9)*lamda(9)+Cm(10)*lamda(10)+Cm(11)*lamda(11)+Cm(12)*lamda(12));
    r_sol = l/A_cell/kappa*100;
    r_hardware = 1; %ohms
    
    %surface area calculation for cathode
    v_cat = Cm(34:41)./rho(2:9);
    vfrac_cat = v_cat/sum(v_cat);
    S_cat_p = (S_cat*vfrac_cat)';
    
    %limiting current calculations
    iLc_cell(1) = z(1)*F*km(1)*Cm(1)+eps;
    iLc_cell(2) = z(2)*F*km(2)*Cm(2)+eps;
    iLc_cell(3) = z(3)*F*km(3)*Cm(3)+eps;
    iLc_cell(4) = z(4)*F*km(4)*Cm(4)+eps;
    iLc_cell(5) = z(5)*F*km(6)*Cm(6)+eps;
    iLc_cell(6) = z(6)*F*km(5)*Cm(5)+eps;
    iLc_cell(7) = z(7)*F*km(7)*Cm(7)+eps;
    iLc_cell(8) = z(8)*F*km(8)*Cm(8)+eps;
    iLc_cell(9) = z(9)*F*km(9)*Cm(9)+eps;
    iLc_cell(10) = z(10)*F*km(10)*Cm(10)+eps;
    iLc_cell(11) = z(11)*F*km(11)*Cm(10)/4+eps;
    
    iLa_cell = -1*ones(1,11);
    iLa_cell(5) = z(5)*F*km(5)*Cm(5)+eps;
    
    foptions = optimoptions(@fsolve, 'Display','off', 'MaxFunctionEvaluations', 3000);
    
    %%%Electrowinning Cell solving%%%
    %solve cell currents and electrode potentials
    onCathode = [1 1 1 1 1 1 1 1 1 1 0];
    onAnode = [0 0 0 0 1 0 0 0 0 0 1];
    if mode == 1
        solver = @(x) cell_solver_p(x(1), x(2), x(3), V_app, r_sol, r_hardware, ...
            Erev_cell, iLa_cell, iLc_cell, onCathode, onAnode, psbl_cell, S_an, S_cat_p, temp);
        %initial guesses [I_an, E_an, E_cat]
        x0 = [0.2, 0.1, -0.1];
        [x, ~, exitflag_cell, ~] = fsolve(solver, x0, foptions);
        I = x(1);
    elseif mode == 2
        solver = @(x) cell_solver_g(x(1), x(2), x(3), I_app, r_sol, r_hardware, ...
            Erev_cell, iLa_cell, iLc_cell, onCathode, onAnode, psbl_cell, S_an, S_cat_p, temp);
        %initial guesses [V, E_an, E_cat]
        x0 = [0.3, 0.1, -0.1];
        [x, ~, exitflag_cell, ~] = fsolve(solver, x0, foptions);
        V = x(1);
    end
    E_an = x(2);
    E_cat = x(3);
    eta_cat = E_cat - Erev_cell;
    eta_an = E_an - Erev_cell;
    i_cat_temp = onCathode.*(i_BV(eta_cat, i0, iLa_cell, iLc_cell, alphas, z, temp));
    i_cat_cat = psbl_cell(1,:).*-subplus(-i_cat_temp);
    i_cat_an = psbl_cell(2,:).*subplus(i_cat_temp);
    I_cat_cat = i_cat_cat*S_cat;
    I_cat_an = i_cat_an.*[S_cat_p(1:4) S_cat S_cat_p(5:8) 0 S_cat];
    I_cat = I_cat_cat+I_cat_an;
    i_an = onAnode.*i_BV(eta_an, i0, iLa_cell, iLc_cell, alphas, z, temp);
    I_an = i_an*S_an;
    I_cell = I_cat+I_an; %overall current for rxn i in cell
    
    %%%Leaching Unit solving%%%
    iLc_corr(1) = z(1)*F*km(1)*Cm(13)+eps;
    iLc_corr(2) = z(2)*F*km(2)*Cm(14)+eps;
    iLc_corr(3) = z(3)*F*km(3)*Cm(15)+eps;
    iLc_corr(4) = z(4)*F*km(4)*Cm(16)+eps;
    iLc_corr(5) = z(5)*F*km(6)*Cm(18)+eps;
    iLc_corr(6) = z(6)*F*km(5)*Cm(17)+eps;
    iLc_corr(7) = z(7)*F*km(7)*Cm(19)+eps;
    iLc_corr(8) = z(8)*F*km(8)*Cm(20)+eps;
    iLc_corr(9) = z(9)*F*km(9)*Cm(21)+eps;
    iLc_corr(10) = z(10)*F*km(10)*Cm(22)+eps;
    iLc_corr(11) = z(11)*F*km(11)*Cm(22)/4+eps;
    
    iLa_corr = -1*ones(1,11);
    iLa_corr(5) = z(5)*F*km(5)*Cm(17)+eps;
    
    %solve extraction lch corrosion rate
    on_PCB_cathode = [1 1 1 1 1 1 1 1 1 1 1];
    on_PCB_anode = [1 1 1 1 1 1 1 1 1 0 1];
    j0 = (Erev_lch(1)+Erev_lch(5))/2; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, Erev_lch, iLa_corr, iLc_corr, S_PCB, ... 
        on_PCB_cathode, on_PCB_anode, psbl_lch, temp);
    %E_corr = fminbnd(cor_solver,min(Erev_lch),max(Erev_lch));
    [E_corr, ~, exitflag_cor, ~] = fsolve(cor_solver, j0, foptions);
    i_BV_corr = i_BV(E_corr-Erev_lch, i0, iLa_corr, iLc_corr, alphas, z, temp);
    i_corr_an = on_PCB_anode.*subplus(i_BV_corr);
    i_corr_cat = psbl_lch.*on_PCB_cathode.*(-subplus(-i_BV_corr));
    i_corr = i_corr_an + i_corr_cat;
    S_corr = [S_PCB(2:5) sum(S_PCB(2:9)) S_PCB(6:9) sum(S_PCB(2:9)) sum(S_PCB(2:9))];
    I_corr = S_corr.*i_corr;
    
    flag(1:11) = E_corr-Erev_lch;
    if exitflag_cell <= 0
        flag(12) = 0;
    end
    if exitflag_cor <= 0
        flag(13) = 0;
    end
end
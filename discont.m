function [flag, isterminal, direction] = discont(t, Cm, temp, pres, vol_cell, ...
    vol_lch, Q, S_an, S_cat, mode, VI_app, n_particles, l, A_cell, solution, foptions)
    direction = [];
    isterminal = zeros(1,13);
    isterminal(12) = 1;
    isterminal(13) = 1;
    flag = ones(1,13);
    disp(['Checking for fsolve failure at t = ' num2str(t)]);
    global F z i0 km alphas lamda rho mw aH2 aO2
    if mode == 1 %potentiostat
        V_app = VI_app;
    elseif mode == 2 %galvanostat
        I_app = VI_app;
    end
    %calculate total mass and wt fractions based on partial masses at
    %timestep
    m_PCB = (Cm(31:37).');
    m_PCB_total = sum(m_PCB);
    %Convert weight to volume and volume fraction
    V_PCB = m_PCB./rho; %Element volumes in m^3
    vfrac_PCB = V_PCB/sum(V_PCB);
    V_PCB_total = sum(V_PCB);
    
    %calculate particle radius in m
    r_particles = (V_PCB_total/(n_particles*4*pi)*3)^(1/3);
    
    %Surface area calculation for corrosion
    SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
    S_PCB_total = V_PCB_total*SSA; %m2 Total surface area of spheres in lch.
    S_PCB = vfrac_PCB*S_PCB_total; %array with exposed surface area of each component
    S_PCB = S_PCB*100^2; %convert m^2 to cm^2
    [Erev_cat, Erev_an, Erev_lch] = nernstPotential(Cm,temp,solution);
    
    %resistance calculation for IR drop in cell
    kappa_cat = 1000*(Cm(1)*lamda(1)+Cm(2)*lamda(2)+Cm(3)*lamda(3)+...
        Cm(4)*lamda(4)+Cm(5)*lamda(5)+Cm(6)*lamda(6)+Cm(7)*lamda(7)+...
        Cm(8)*lamda(8)+Cm(9)*lamda(9)+Cm(10)*lamda(10));
    kappa_an = 1000*(Cm(11)*lamda(1)+Cm(12)*lamda(2)+Cm(13)*lamda(3)+...
        Cm(14)*lamda(4)+Cm(15)*lamda(5)+Cm(16)*lamda(6)+Cm(17)*lamda(7)+...
        Cm(18)*lamda(8)+Cm(19)*lamda(9)+Cm(20)*lamda(10));
    kappa_avg = mean([kappa_an kappa_cat]);
    r_sol = l/A_cell/kappa_avg*100;
    r_hardware = 10; %ohms
    
    %cell volume division
    vol_cat = vol_cell/2;
    vol_an = vol_cat;
    %surface area calculation for cathode
    v_cat = Cm(38:43)./rho(2:7);
    vfrac_cat = v_cat/sum(v_cat);
    S_cat_p = (S_cat*vfrac_cat)';
    
    %limiting current calculations for cathode side
    iLc_cat(1) = z(1)*F*km(1)*Cm(1)+eps;
    iLc_cat(2) = z(2)*F*km(2)*Cm(2)+eps;
    iLc_cat(3) = z(3)*F*km(4)*Cm(4)+eps;
    iLc_cat(4) = z(4)*F*km(3)*Cm(3)+eps;
    iLc_cat(6) = -1;
    iLc_cat(8) = z(8)*F*km(10)*Cm(10)+eps;
    iLc_cat(10) = z(10)*F*km(8)*Cm(8)+eps;
    iLc_cat(11) = z(11)*F*km(8)*Cm(8)+eps;
    
    iLa_cat = -1*ones(1,11);
    iLa_cat(3) = z(5)*F*km(3)*Cm(3)+eps;
    iLa_cat(6) = z(6)*F*km(9)*Cm(9)+eps;
    iLa_cat(8) = z(8)*F*km(9)*Cm(9)/4+eps;
    
    if solution == 1 %Cl-, base metal system
    	iLc_cat(5) = z(5)*F*km(5)*Cm(5)+eps;
        iLc_cat(7) = z(7)*F*km(6)*Cm(6)+eps;
        iLc_cat(9) = z(9)*F*km(7)*Cm(7)+eps;
    else %S2O3, precious metal system
        iLc_cat(5) = z(5)*F*km(5)*Cm(5)+eps;
        iLc_cat(7) = z(7)*F*km(6)*Cm(6)+eps;
        iLc_cat(9) = z(9)*F*km(7)*Cm(7)+eps;
        iLa_cat(5) = z(5)*F*km(9)*Cm(9)/2+eps;
        iLa_cat(7) = z(7)*F*km(9)*Cm(9)/2+eps;
        iLa_cat(9) = z(9)*F*km(9)*Cm(9)/4+eps;
    end
    
    %limiting current calculations for anode side
    iLc_an(1) = z(1)*F*km(1)*Cm(11)+eps;
    iLc_an(2) = z(2)*F*km(2)*Cm(12)+eps;
    iLc_an(3) = z(3)*F*km(4)*Cm(14)+eps;
    iLc_an(4) = z(4)*F*km(3)*Cm(13)+eps;
    iLc_an(6) = -1;
    iLc_an(8) = z(8)*F*km(10)*Cm(20)+eps;
    iLc_an(10) = z(10)*F*km(8)*Cm(18)+eps;
    iLc_an(11) = z(11)*F*km(8)*Cm(18)+eps;
    
    iLa_an = -1*ones(1,11);
    iLa_an(3) = z(5)*F*km(3)*Cm(13)+eps;
    iLa_an(6) = z(6)*F*km(9)*Cm(19)+eps;
    iLa_an(8) = z(8)*F*km(9)*Cm(19)/4+eps;
    
    if solution == 1 %Cl-, base metal system
    	iLc_an(5) = z(5)*F*km(5)*Cm(15)+eps;
        iLc_an(7) = z(7)*F*km(6)*Cm(16)+eps;
        iLc_an(9) = z(9)*F*km(7)*Cm(17)+eps;
    else %S2O3, precious metal system
        iLc_an(5) = z(5)*F*km(5)*Cm(15)+eps;
        iLc_an(7) = z(7)*F*km(6)*Cm(16)+eps;
        iLc_an(9) = z(9)*F*km(7)*Cm(17)+eps;
        iLa_an(5) = z(5)*F*km(9)*Cm(19)/2+eps;
        iLa_an(7) = z(7)*F*km(9)*Cm(19)/2+eps;
        iLa_an(9) = z(9)*F*km(9)*Cm(19)/4+eps;
    end
    
    %%%Electrowinning Cell solving%%%
    %solve cell currents and electrode potentials
    onCathode = [1 1 1 1 1 0 1 1 1 1 1];
    onAnode = [0 0 1 0 0 0 0 0 0 0 1];
    if mode == 1
        solver = @(x) cell_solver_p(x(1), x(2), x(3), V_app, r_sol, r_hardware, ...
            Erev_cat, Erev_an, iLa_cat, iLc_cat, iLa_an, iLa_cat, onCathode, ...
            onAnode, S_an, S_cat_p, temp);
        %initial guesses [I_an, E_an, E_cat]
        x0 = [0.2, 0.1, -0.1];
        [x,~,exitflag_cell,o_cell] = fsolve(solver, x0, foptions);
        I = x(1);
    elseif mode == 2
        solver = @(x) cell_solver_g(x(1), x(2), x(3), I_app, r_sol, r_hardware, ...
            Erev_cat, Erev_an, iLa_cat, iLc_cat, iLa_an, iLa_cat, onCathode, ...
            onAnode, S_an, S_cat_p, temp);
        %initial guesses [V, E_an, E_cat]
        x0 = [2, 0.5, -0.5];
        [x,~,exitflag_cell,o_cell] = fsolve(solver, x0, foptions);
        V = x(1);
    end
    E_an = x(2);
    E_cat = x(3);
    eta_cat = E_cat - Erev_cat;
    eta_an = E_an - Erev_an;
    i_cat = onCathode.*(i_BV(eta_cat, i0, iLa_cat, iLc_cat, alphas, z, temp));
    i_cat_cat = -subplus(-i_cat);
    i_cat_an = subplus(i_cat);
    I_cat_cat = i_cat_cat*S_cat;
    I_cat_cat(6) = 0; %silver can't be plated from AgCl(s)
    I_cat_an = i_cat_an.*[S_cat_p(1:2) S_cat S_cat_p(3:4) S_cat_p(4:5) S_cat_p(5:6) 0 S_cat];
    I_cat = I_cat_cat+I_cat_an;
    i_an = onAnode.*i_BV(eta_an, i0, iLa_an, iLc_an, alphas, z, temp);
    I_an = i_an*S_an;
    I_cell = I_cat+I_an; %overall current for rxn i in cell
    
    %%%Leaching Unit solving%%%
    iLc_corr(1) = z(1)*F*km(1)*Cm(21)+eps;
    iLc_corr(2) = z(2)*F*km(2)*Cm(22)+eps;
    iLc_corr(3) = z(3)*F*km(4)*Cm(24)+eps;
    iLc_corr(4) = z(4)*F*km(3)*Cm(23)+eps;
    iLc_corr(6) = -1;
    iLc_corr(8) = z(8)*F*km(10)*Cm(30)+eps;
    iLc_corr(10) = z(10)*F*km(8)*Cm(18)+eps;
    iLc_corr(11) = z(11)*F*km(8)*Cm(18)+eps;
    
    iLa_corr = -1*ones(1,11);
    iLa_corr(3) = z(3)*F*km(3)*Cm(23)+eps;
    
    if solution == 1
        iLc_corr(5) = z(5)*F*km(5)*Cm(25)+eps;
        iLc_corr(7) = z(7)*F*km(6)*Cm(26)+eps;
        iLc_corr(9) = z(9)*F*km(7)*Cm(27)+eps;
    else
        iLc_corr(5) = z(5)*F*km(5)*Cm(25)+eps;
        iLc_corr(7) = z(7)*F*km(6)*Cm(26)+eps;
        iLc_corr(9) = z(9)*F*km(7)*Cm(27)+eps;
        iLa_corr(5) = z(5)*F*km(9)*Cm(29)/2+eps;
        iLa_corr(7) = z(7)*F*km(9)*Cm(29)/2+eps;
        iLa_corr(9) = z(9)*F*km(9)*Cm(29)/4+eps;
    end
    
    %solve extraction lch corrosion rate
    on_PCB_cathode = [1 1 1 1 1 0 1 1 1 1 1];
    on_PCB_anode = [1 1 1 1 1 1 1 1 1 0 1];
    j0 = (Erev_lch(1)+Erev_lch(5))/2; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, Erev_lch, iLa_corr, iLc_corr, S_PCB, ... 
        on_PCB_cathode, on_PCB_anode, temp);
    %E_corr = fminbnd(cor_solver,min(Erev_lch),max(Erev_lch));
    [E_corr,~,exitflag_cor,o_cor] = fsolve(cor_solver, j0, foptions);
    i_BV_corr = i_BV(E_corr-Erev_lch, i0, iLa_corr, iLc_corr, alphas, z, temp);
    i_corr_an = on_PCB_anode.*subplus(i_BV_corr);
    i_corr_cat = on_PCB_cathode.*(-subplus(-i_BV_corr));
    i_corr = i_corr_an + i_corr_cat;
    %arrange surface areas in proper order
    S_corr = [S_PCB(2:3) sum(S_PCB(2:7)) S_PCB(4:5) S_PCB(5:6) S_PCB(6:7) sum(S_PCB(2:7)) sum(S_PCB(2:7))];
    I_corr = S_corr.*i_corr;
    
    flag(1:11) = E_corr-Erev_lch;
    if exitflag_cell <= 0
        flag(12) = 0;
        disp(['Electrowinning cell solver failed. Message: ' o_cell.message]);
    end
    if exitflag_cor <= 0
        flag(13) = 0;
        disp(['Leaching solver failed. Message: ' o_cor.message]);
    end
end
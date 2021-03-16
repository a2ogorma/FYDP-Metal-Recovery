function [flag, isterminal, direction] = discont(t, Cm, temp, pres, vol_cell, ...
    vol_lch, Q, S_an, S_cat, mode, VI_app, n_particles, l, A_cell, n_units, solution, iL_default, foptions)
    direction = [];
    isterminal = zeros(1,13);
    isterminal(12) = 0;
    isterminal(13) = 0;
    flag = ones(1,13);
    disp(['Checking for fsolve failure at t = ' num2str(t)]);
    global F z km_cell lamda rho rho_e mu_e Dab Sc 
    
    if mode == 1 %potentiostat
        V_app = VI_app;
    elseif mode == 2 %galvanostat
        I_app = VI_app/n_units;
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
    vol_unit = vol_cell/n_units;
    vol_cat = vol_unit/2;
    vol_an = vol_cat;
    %surface area calculation for cathode
    v_cat = Cm(38:43)./rho(2:7);
    vfrac_cat = v_cat/sum(v_cat);
    S_cat_p = (S_cat*vfrac_cat)';
    
    %%%Electrowinning Cell solving%%%
    %limiting current calculations for cathode side
    iLc_cat(1) = z(1)*F*km_cell(1)*Cm(1)+eps;
    iLc_cat(2) = z(2)*F*km_cell(2)*Cm(2)+eps;
    iLc_cat(3) = z(3)*F*km_cell(4)*Cm(4)+eps;
    iLc_cat(4) = z(4)*F*km_cell(3)*Cm(3)+eps;
    iLc_cat(6) = iL_default;
    iLc_cat(8) = z(8)*F*km_cell(10)*Cm(10)+eps;
    iLc_cat(10) = z(10)*F*km_cell(8)*Cm(8)+eps;
    iLc_cat(11) = z(11)*F*km_cell(8)*Cm(8)+eps;
    
    iLa_cat = iL_default*ones(1,11);
    iLa_cat(3) = z(5)*F*km_cell(3)*Cm(3)+eps;
    iLa_cat(6) = z(6)*F*km_cell(9)*Cm(9)+eps;
    iLa_cat(8) = z(8)*F*km_cell(9)*Cm(9)/4+eps;
    
    if solution == 1 %Cl-, base metal system
    	iLc_cat(5) = z(5)*F*km_cell(5)*Cm(5)+eps;
        iLc_cat(7) = z(7)*F*km_cell(6)*Cm(6)+eps;
        iLc_cat(9) = z(9)*F*km_cell(7)*Cm(7)+eps;
    else %S2O3, precious metal system
        iLc_cat(5) = z(5)*F*km_cell(5)*Cm(5)+eps;
        iLc_cat(7) = z(7)*F*km_cell(6)*Cm(6)+eps;
        iLc_cat(9) = z(9)*F*km_cell(7)*Cm(7)+eps;
        iLa_cat(5) = z(5)*F*km_cell(9)*Cm(9)/2+eps;
        iLa_cat(7) = z(7)*F*km_cell(9)*Cm(9)/2+eps;
        iLa_cat(9) = z(9)*F*km_cell(9)*Cm(9)/4+eps;
    end
    
    %limiting current calculations for anode side
    iLc_an(1) = z(1)*F*km_cell(1)*Cm(11)+eps;
    iLc_an(2) = z(2)*F*km_cell(2)*Cm(12)+eps;
    iLc_an(3) = z(3)*F*km_cell(4)*Cm(14)+eps;
    iLc_an(4) = z(4)*F*km_cell(3)*Cm(13)+eps;
    iLc_an(6) = iL_default;
    iLc_an(8) = z(8)*F*km_cell(10)*Cm(20)+eps;
    iLc_an(10) = z(10)*F*km_cell(8)*Cm(18)+eps;
    iLc_an(11) = z(11)*F*km_cell(8)*Cm(18)+eps;
    
    iLa_an = iL_default*ones(1,11);
    iLa_an(3) = z(5)*F*km_cell(3)*Cm(13)+eps;
    iLa_an(6) = z(6)*F*km_cell(9)*Cm(19)+eps;
    iLa_an(8) = z(8)*F*km_cell(9)*Cm(19)/4+eps;
    
    if solution == 1 %Cl-, base metal system
    	iLc_an(5) = z(5)*F*km_cell(5)*Cm(15)+eps;
        iLc_an(7) = z(7)*F*km_cell(6)*Cm(16)+eps;
        iLc_an(9) = z(9)*F*km_cell(7)*Cm(17)+eps;
    else %S2O3, precious metal system
        iLc_an(5) = z(5)*F*km_cell(5)*Cm(15)+eps;
        iLc_an(7) = z(7)*F*km_cell(6)*Cm(16)+eps;
        iLc_an(9) = z(9)*F*km_cell(7)*Cm(17)+eps;
        iLa_an(5) = z(5)*F*km_cell(9)*Cm(19)/2+eps;
        iLa_an(7) = z(7)*F*km_cell(9)*Cm(19)/2+eps;
        iLa_an(9) = z(9)*F*km_cell(9)*Cm(19)/4+eps;
    end
    %Convert units to A/cm^2 From A*m/dm^3 
    iLc_an = iLc_an*0.1;
    iLc_cat = iLc_cat*0.1;
    iLa_an = iLa_an*0.1;
    iLa_cat = iLa_cat*0.1;
    
    %solve cell currents and electrode potentials
    onCathode = [1 1 1 1 1 0 1 1 1 1 1];
    onAnode = [0 0 1 0 0 0 0 0 0 0 1];
    global x0 %initial guesses [I_an, E_an, E_cat] for mode 1 or [V, E_an, E_cat] for mode 2
    if mode == 1
        solver = @(x) cell_solver_p(x(1), x(2), x(3), V_app, r_sol, r_hardware, ...
            Erev_cat, Erev_an, iLa_cat, iLc_cat, iLa_an, iLc_an, onCathode, ...
            onAnode, S_an, S_cat_p, temp);
        [x,~,exitflag_cell,o_cell] = fsolve(solver, x0, foptions);
        I = x(1);
    elseif mode == 2
        solver = @(x) cell_solver_g(x(1), x(2), x(3), I_app, r_sol, r_hardware, ...
            Erev_cat, Erev_an, iLa_cat, iLc_cat, iLa_an, iLc_an, onCathode, ...
            onAnode, S_an, S_cat_p, temp);
        [x,~,exitflag_cell,o_cell] = fsolve(solver, x0, foptions);
        V = x(1);
    end
    
    %%%Leaching Unit solving%%%
    u_lch = 0.5; %m/s assumed in stirred tank
    Re_lch = rho_e*u_lch*r_particles*2/mu_e;
    Pe_lch = Re_lch.*Sc;
    Sh_lch = (4+1.21*Pe_lch.^(2/3)).^0.5;
    km_lch = Dab.*Sh_lch/r_particles/2;
    
    iLc_corr(1) = z(1)*F*km_lch(1)*Cm(21)+eps;
    iLc_corr(2) = z(2)*F*km_lch(2)*Cm(22)+eps;
    iLc_corr(3) = z(3)*F*km_lch(4)*Cm(24)+eps;
    iLc_corr(4) = z(4)*F*km_lch(3)*Cm(23)+eps;
    iLc_corr(6) = iL_default;
    iLc_corr(8) = z(8)*F*km_lch(10)*Cm(30)+eps;
    iLc_corr(10) = z(10)*F*km_lch(8)*Cm(18)+eps;
    iLc_corr(11) = z(11)*F*km_lch(8)*Cm(18)+eps;
    
    iLa_corr = iL_default*ones(1,11);
    iLa_corr(3) = z(3)*F*km_lch(3)*Cm(23)+eps;
    
    if solution == 1
        iLc_corr(5) = z(5)*F*km_lch(5)*Cm(25)+eps;
        iLc_corr(7) = z(7)*F*km_lch(6)*Cm(26)+eps;
        iLc_corr(9) = z(9)*F*km_lch(7)*Cm(27)+eps;
    else
        iLc_corr(5) = z(5)*F*km_lch(5)*Cm(25)+eps;
        iLc_corr(7) = z(7)*F*km_lch(6)*Cm(26)+eps;
        iLc_corr(9) = z(9)*F*km_lch(7)*Cm(27)+eps;
        iLa_corr(5) = z(5)*F*km_lch(9)*Cm(29)/2+eps;
        iLa_corr(7) = z(7)*F*km_lch(9)*Cm(29)/2+eps;
        iLa_corr(9) = z(9)*F*km_lch(9)*Cm(29)/4+eps;
    end
    %Convert units to A/cm^2 from A*m/dm^3
    iLc_corr = iLc_corr*0.1;
    iLa_corr = iLa_corr*0.1;
    
    %solve extraction lch corrosion rate
    on_PCB_cathode = [1 1 1 1 1 0 1 1 1 1 1];
    on_PCB_anode = [1 1 1 1 1 1 1 1 1 0 1];
    cor_solver = @(E_corr)cor(E_corr, Erev_lch, iLa_corr, iLc_corr, S_PCB, ... 
        on_PCB_cathode, on_PCB_anode, temp);
    global E_corr0
    [E_corr,~,exitflag_cor,o_cor] = fsolve(cor_solver, E_corr0, foptions);
    %{
    j0 = (randperm(20)-1)/19*(max(Erev_lch)-min(Erev_lch))+min(Erev_lch); %Initial guess for E_corr, V b/w max and min Nernst potentials
    for k = 1:1:numel(j0)
        [E_corr,~,exitflag_cor,o_cor] = fsolve(cor_solver, j0(k), foptions);
        if exitflag_cor >= 1
            break
        end
    end
    %}
    E_corr0 = E_corr; %Set new initial guess to old one
    flag(1:11) = E_corr-Erev_lch;
    if exitflag_cell <= 0
        flag(12) = 0;
        disp(['Electrowinning cell solver failed. Message: ' o_cell.message]);
    end
    if exitflag_cor <= 0
        flag(13) = 0;
        disp(['Leaching solver failed. Message: ' o_cor.message 'Cor result:' num2str(cor_solver(E_corr))]);
        %{
        e_corr = [E_corr-0.5:0.01:E_corr+0.5];
        for i = 1:1:numel(e_corr)
            c(i) = cor_solver(e_corr(i));
        end
        plot(e_corr, c,'b-', E_corr, cor_solver(E_corr), 'r.');
        %}
    end
end
function dt = ion_balance(t, Cm, temp, pres, vol_cell, vol_lch, Q, S_an, ...
    S_cat, mode, VI_app, n_particles, l, A_cell, solution, iL_default, foptions)
    %t = time span (s)
    %Cm = ionic concentration and solid mass vector (M, kg)
    %pres = system pressure (atm)
    %vol_cell = electrowinning cell volume (L)
    %vol_lch = leaching vessel volume (L)
    %Q = circulation flowrate (L/s)
    %S_an = cell anodic surface area (cm^2);
    %S_cat = cell cathodic surface area (cm^2);
    %V_app = applied voltage to cell (V);
    %n_particles = number of PCB particles as calculated by radius/volume;
    %l = cell distance between electrodes;
    %A_cell = cell cross sectional area;
    %{
    Ionic species order:
    Base: Cu2+, Sn2+, Fe2+, Fe3+, Ag+, Au3+, Pd2+, H+, Cl-, (AuCl4)-,
    Precious: Cu2+, Sn2+, Fe2+, Fe3+, [Ag(S2O3)2]3-, [Au(S2O3)2]3-,
        [Pd(S2O3)4]6-, H+, (S2O3)2-, (AuCl4)-
    Solid mass order: Inert (If applicable), Cu, Sn, Fe, Ag, Au, Pd

    Reactions
    Cu2+ + 2e- <--> Cu(s) (1)
    Sn2+ + 2e- <--> Sn(s) (2)
    Fe3+ + e- <--> Fe2+ (3)
    Fe2+ + 2e- <--> Fe(s) (4)
    Ag+ + e- <--> Ag(s) (5B) OR [Ag(S2O3)2]3- + e- <--> Ag(s) + 2(S2O3)2- (5P)
    AgCl(s) + e- <--> Ag(s) + Cl- (6)
    Au3+ + 3e- <--> Au(s) (7B) OR [Au(S2O3)2]3- + e- <--> Au(s) + 2(S2O3)2- (7P)
    [AuCl4]- + 3e- <--> Au(s) + 4Cl- (8)
    Pd2+ + 2e- <--> Pd(s) (9B) OR [Pd(S2O3)4]6- + 2e- <--> Pd(s) + 4(S2O3)2- (9P)
    2H+ + e- <--> H2(g) (10)
    4H+ + O2(g) + 4e- <--> 2H2O(l) (11)
    %}
    
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
    iLc_cat(6) = iL_default;
    iLc_cat(8) = z(8)*F*km(10)*Cm(10)+eps;
    iLc_cat(10) = z(10)*F*km(8)*Cm(8)+eps;
    iLc_cat(11) = z(11)*F*km(8)*Cm(8)+eps;
    
    iLa_cat = iL_default*ones(1,11);
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
    iLc_an(6) = iL_default;
    iLc_an(8) = z(8)*F*km(10)*Cm(20)+eps;
    iLc_an(10) = z(10)*F*km(8)*Cm(18)+eps;
    iLc_an(11) = z(11)*F*km(8)*Cm(18)+eps;
    
    iLa_an = iL_default*ones(1,11);
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
    global x0 %initial guesses [I_an, E_an, E_cat] for mode 1 or [V, E_an, E_cat] for mode 2
    if mode == 1
        solver = @(x) cell_solver_p(x(1), x(2), x(3), V_app, r_sol, r_hardware, ...
            Erev_cat, Erev_an, iLa_cat, iLc_cat, iLa_an, iLc_an, onCathode, ...
            onAnode, S_an, S_cat_p, temp);
        x = fsolve(solver, x0, foptions);
    elseif mode == 2
        solver = @(x) cell_solver_g(x(1), x(2), x(3), I_app, r_sol, r_hardware, ...
            Erev_cat, Erev_an, iLa_cat, iLc_cat, iLa_an, iLc_an, onCathode, ...
            onAnode, S_an, S_cat_p, temp);
        %initial guesses 
        x = fsolve(solver, x0, foptions);
    end
    x0 = x; %set initial guess for next step to computed solution
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
    iLc_corr(6) = iL_default;
    iLc_corr(8) = z(8)*F*km(10)*Cm(30)+eps;
    iLc_corr(10) = z(10)*F*km(8)*Cm(18)+eps;
    iLc_corr(11) = z(11)*F*km(8)*Cm(18)+eps;
    
    iLa_corr = iL_default*ones(1,11);
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
    cor_solver = @(E_corr)cor(E_corr, Erev_lch, iLa_corr, iLc_corr, S_PCB, ... 
        on_PCB_cathode, on_PCB_anode, temp);
    j0 = (randperm(20)-1)/19*(max(Erev_lch)-min(Erev_lch))+min(Erev_lch); %Initial guess for E_corr, V b/w max and min Nernst potentials
    for k = 1:1:numel(j0)
        [E_corr,~,exitflag_cor,~] = fsolve(cor_solver, j0(k), foptions);
        if exitflag_cor >= 1
            break
        end
    end
    %E_corr = fminbnd(cor_solver,min(Erev_lch),max(Erev_lch));
    i_BV_corr = i_BV(E_corr-Erev_lch, i0, iLa_corr, iLc_corr, alphas, z, temp);
    i_corr_an = on_PCB_anode.*subplus(i_BV_corr);
    i_corr_cat = on_PCB_cathode.*(-subplus(-i_BV_corr));
    i_corr = i_corr_an + i_corr_cat;
    %arrange surface areas in proper order
    S_corr = [S_PCB(2:3) sum(S_PCB(2:7)) S_PCB(4:5) S_PCB(5:6) S_PCB(6:7) sum(S_PCB(2:7)) sum(S_PCB(2:7))];
    I_corr = S_corr.*i_corr;
    
    %%%Solving for mass balance differentials%%%%
    %Calculate concentration/mass balances
    dt = zeros(size(Cm));
    %Electrowinning cathodic side concentration balances
    dt(1) = ((Cm(21)-Cm(1))*Q + I_cat(1)/F/z(1))/vol_cat; %Cu2+
    dt(2) = ((Cm(22)-Cm(2))*Q + I_cat(2)/F/z(2))/vol_cat; %Sn2+
    dt(3) = ((Cm(23)-Cm(3))*Q + (I_cat(4)/F/z(4))-(I_cat(3)/F/z(3)))/vol_cat; %Fe2+
    dt(4) = ((Cm(24)-Cm(4))*Q + I_cat(3)/F/z(3))/vol_cat; %Fe3+
    dt(5) = ((Cm(25)-Cm(5))*Q + I_cat(5)/F/z(5))/vol_cat; %Ag+ or (AgS2O3)3-
    dt(6) = ((Cm(26)-Cm(6))*Q + I_cat(7)/F/z(7))/vol_cat; %Au3+ or (AuS2O3)3-
    dt(7) = ((Cm(27)-Cm(7))*Q + I_cat(9)/F/z(9))/vol_cat; %Pd2+ or (PdS2O3)6-
    dt(8) = ((Cm(28)-Cm(8))*Q + I_cat(10)/F/z(10)+I_cat(11)/F/z(11))/vol_cat; %H+
    if solution == 1
        dt(9) = ((Cm(29)-Cm(9))*Q - I_cat(6)/F/z(6)-I_cat(8)*4/F/z(8))/vol_cat; %Cl-
    else
        dt(9) = ((Cm(29)-Cm(9))*Q - I_cat(5)*2/F/z(5) - I_cat(7)*2/F/z(7)...
            - I_cat(9)*4/F/z(9))/vol_cat; %(S203)2-
    end
    dt(10) = ((Cm(30)-Cm(10))*Q + I_cat(10)/F/z(10) + I_cat(11)/F/z(11))/vol_cat; %AuCl4-
   
    %Electrowinning anodic side concentration balances
    dt(11) = ((Cm(1)-Cm(11))*Q + I_an(1)/F/z(1))/vol_an; %Cu2+
    dt(12) = ((Cm(2)-Cm(12))*Q + I_an(2)/F/z(2))/vol_an; %Sn2+
    dt(13) = ((Cm(3)-Cm(13))*Q + I_an(4)/F/z(4)-(I_an(3)/F/z(3)))/vol_an; %Fe2+
    dt(14) = ((Cm(4)-Cm(14))*Q + I_an(3)/F/z(3))/vol_an; %Fe3+
    dt(15) = ((Cm(5)-Cm(15))*Q + I_an(5)/F/z(5))/vol_an; %Ag+ or (AgS2O3)3-
    dt(16) = ((Cm(6)-Cm(16))*Q + I_an(7)/F/z(7))/vol_an; %Au3+ or (AuS2O3)3-
    dt(17) = ((Cm(7)-Cm(17))*Q + I_an(9)/F/z(9))/vol_an; %Pd2+ or (PdS2O3)6-
    dt(18) = ((Cm(8)-Cm(18))*Q + I_an(10)/F/z(10)+I_an(11)/F/z(11))/vol_an; %H+
    if solution == 1
        dt(19) = ((Cm(9)-Cm(19))*Q - I_an(6)/F/z(6)-I_an(8)*4/F/z(8))/vol_an; %Cl-
    else
        dt(19) = ((Cm(9)-Cm(19))*Q - I_an(5)*2/F/z(5) - I_an(7)*2/F/z(7)...
            - I_an(9)*4/F/z(9))/vol_an; %(S203)2-
    end
    dt(20) = ((Cm(10)-Cm(20))*Q + I_an(10)/F/z(10) + I_an(11)/F/z(11))/vol_an; %AuCl4-
    
    %Leaching vessel concentration balances
    dt(21) = ((Cm(11)-Cm(21))*Q + I_corr(1)/F/z(1))/vol_lch; %Cu2+
    dt(22) = ((Cm(12)-Cm(22))*Q + I_corr(2)/F/z(2))/vol_lch; %Sn2+
    dt(23) = ((Cm(13)-Cm(23))*Q + I_corr(4)/F/z(4) - I_corr(3)/F/z(3))/vol_lch; %Fe2+
    dt(24) = ((Cm(14)-Cm(24))*Q + I_corr(3)/F/z(3))/vol_lch; %Fe3+
    dt(25) = ((Cm(15)-Cm(25))*Q + I_corr(5)/F/z(5))/vol_lch; %Ag+ or (AgS2O3)3-
    dt(26) = ((Cm(16)-Cm(26))*Q + I_corr(7)/F/z(7))/vol_lch; %Au3+ or (AuS2O3)3-
    dt(27) = ((Cm(17)-Cm(27))*Q + I_corr(9)/F/z(9))/vol_lch; %Pd2+ or (PdS2O3)6-
    dt(28) = ((Cm(18)-Cm(28))*Q + I_corr(10)/F/z(10) + I_corr(11)/F/z(11))/vol_lch; %H+
    if solution == 1
        dt(29) = ((Cm(19)-Cm(29))*Q - I_corr(6)/F/z(6)-I_corr(8)*4/F/z(8))/vol_an; %Cl-
    else
        dt(29) = ((Cm(19)-Cm(29))*Q - I_corr(5)*2/F/z(5) - I_corr(7)*2/F/z(7)...
            - I_corr(9)*4/F/z(9))/vol_an; %(S203)2-
    end
    dt(30) = ((Cm(20)-Cm(30))*Q + I_corr(10)/F/z(10) + I_corr(11)/F/z(11))/vol_an; %AuCl4-
    
    %PCB metal mass balances in kg units
    dt(31) = 0; %Inert material 
    dt(32) = -mw(2)*I_corr(1)/F/z(1)/1000; %Cu(s)
    dt(33) = -mw(3)*I_corr(2)/F/z(2)/1000; %Sn(s)
    dt(34) = -mw(4)*I_corr(4)/F/z(4)/1000; %Fe(s)
    dt(35) = -mw(5)*(I_corr(5)/F/z(5)+I_corr(6)/F/z(6))/1000; %Ag(s)
    dt(36) = -mw(6)*(I_corr(7)/F/z(7)+I_corr(8)/F/z(8))/1000; %Au(s)
    dt(37) = -mw(7)*I_corr(9)/F/z(9)/1000; %Pd(s)
    
    %Recovered metals on cathode
    dt(38) = -mw(2)*I_cat(1)/F/z(1)/1000; %Cu(s)
    dt(39) = -mw(3)*I_cat(2)/F/z(2)/1000; %Sn(s)
    dt(40) = -mw(4)*I_cat(4)/F/z(4)/1000; %Fe(s)
    dt(41) = -mw(5)*(I_cat(5)/F/z(5)+I_cat(6)/F/z(6))/1000; %Ag(s)
    dt(42) = -mw(6)*(I_cat(7)/F/z(7)+I_cat(8)/F/z(8))/1000; %Au(s)
    dt(43) = -mw(7)*I_cat(9)/F/z(9)/1000; %Pd(s)
    %display current time step
    %t
end
function dt = ion_balance(t, Cm, temp, pres, vol_cell, vol_lch, Q, S_an, S_cat, V_app, n_particles, l, A_cell)
    %t = time span (s)
    %Cm = ionic concentration and solid mass vector (M, kg)
    %Concentration order: Cu2 Sn2 Al3 Pb2 Fe2 Fe3 Ag Au Pd2 H Cl
    %Solid mass order: Inert Cu Sn Al Pb Fe Ag Au Pd
    %pres = system pressure (atm)
    %vol_cell = electrowinning cell volume (L)
    %vol_lch = leaching vessel volume (L)
    %Q = circulation flowrate (L/s)
    %S_an = cell anodic surface area (cm^2);
    %Sn_cat = cell cathodic surface area (cm^2);
    %V_app = applied voltage to cell (V);
    %n_particles = number of PCB particles as calculated by radius/volume;
    %l = cell distance between electrodes;
    %A_cell = cell cross sectional area;
    global R F e z i0 iL km alphas Eo gamma lamda rho mw aH2 aO2
    %{
    Reactions
    Cu2+ + 2e- <--> Cu(s) (1)
    Sn2+ + 2e- <--> Sn(s) (2)
    Al3+ + 3e- <--> Al(s) (3)
    Pb2+ + 2e- <--> Pb(s) (4)
    Fe3+ + e- <--> Fe2+ (5) (Fe1)
    Fe2+ + 2e- <--> Fe(s) (6) (Fe2)
    Ag+ + e- <--> Ag(s) (7)
    Au+ + e- <--> Au(s) (8)
    Pd2+ + 2e- <--> Pd(s) (9)
    2H+ + e- <--> H2(g) (10)
    4H+ + O2(g) + 4e- <--> 2H2O(l) (11)
    %}
    %calculate total mass and wt fractions based on partial masses at
    %timestep
    m_PCB = (Cm(25:33).');
    m_PCB_total = sum(m_PCB);
    wtfrac_PCB = m_PCB/m_PCB_total;
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
    
    [Erev_cell, Erev_lch] = nernstPotential(Cm,temp,aH2,aO2);
    
    %resistance calculation for IR drop in cell
    kappa = 1000*(Cm(1)*lamda(1)+Cm(2)*lamda(2)+Cm(3)*lamda(3)+...
        Cm(4)*lamda(4)+Cm(5)*lamda(5)+Cm(6)*lamda(6)+Cm(7)*lamda(7)+...
        Cm(8)*lamda(8)+Cm(9)*lamda(9)+Cm(10)*lamda(10)+Cm(11)*lamda(11)+Cm(12)*lamda(12));
    r_sol = l/A_cell/kappa*100;
    r_hardware = 1; %ohms
    
    %limiting current calculations
    iL_cell(1) = z(1)*F*km(1)*Cm(1)+eps;
    iL_cell(2) = z(2)*F*km(2)*Cm(2)+eps;
    iL_cell(3) = z(3)*F*km(3)*Cm(3)+eps;
    iL_cell(4) = z(4)*F*km(4)*Cm(4)+eps;
    iL_cell(5) = z(5)*F*km(5)*Cm(6)+eps;
    iL_cell(6) = z(6)*F*km(6)*Cm(5)+eps;
    iL_cell(7) = z(7)*F*km(7)*Cm(7)+eps;
    iL_cell(8) = z(8)*F*km(8)*Cm(8)+eps;
    iL_cell(9) = z(9)*F*km(9)*Cm(9)+eps;
    iL_cell(10) = z(10)*F*km(10)*Cm(10)+eps;
    iL_cell(11) = z(11)*F*km(11)*Cm(10)/4+eps;
    
    %%%Electrowinning Cell solving%%%
    %solve cell currents and electrode potentials
    onCathode = [1 1 1 1 1 1 1 1 1 1 0];
    onAnode = [0 0 0 0 1 0 0 0 0 0 1];
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol, r_hardware,...
        Erev_cell, iL_cell, onCathode, onAnode, S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [0.2, 0.1, -0.1];
    options = optimoptions(@fsolve, 'Display','final', 'MaxFunctionEvaluations', 3000);
    [x, fval] = fsolve(solver, x0, options);
    I = x(1);
    E_an = x(2);
    E_cat = x(3);
    eta_cat = E_cat - Erev_cell;
    eta_an = E_an - Erev_cell;
    i_cat = onCathode.*i_BV(eta_cat, i0, iL_cell, alphas, z, temp);
    I_cat = i_cat*S_cat;
    i_an = onAnode.*i_BV(eta_an, i0, iL_cell, alphas, z, temp);
    I_an = i_an*S_an;
    I_cell = I_cat+I_an; %overall current for rxn i in cell
    
    %%%Leaching Unit solving%%%
    iL_corr(1) = z(1)*F*km(1)*Cm(23)+eps;
    iL_corr(2) = z(2)*F*km(2)*Cm(23)+eps;
    iL_corr(3) = z(3)*F*km(3)*Cm(23)+eps;
    iL_corr(4) = z(4)*F*km(4)*Cm(23)+eps;
    iL_corr(5) = z(5)*F*km(5)*Cm(17)+eps;
    iL_corr(6) = z(6)*F*km(6)*Cm(18)+eps;
    iL_corr(7) = z(7)*F*km(7)*Cm(24)+eps;
    iL_corr(8) = z(8)*F*km(8)*Cm(24)+eps;
    iL_corr(9) = z(9)*F*km(9)*Cm(24)+eps;
    iL_corr(10) = z(10)*F*km(10)*Cm(22)+eps;
    iL_corr(11) = z(11)*F*km(11)*Cm(22)/4+eps;
    
    %solve extraction lch corrosion rate
    on_PCB_cathode = [1 1 1 1 1 1 1 1 1 1 1];
    on_PCB_anode = [1 1 1 1 1 1 1 1 1 0 1];
    fz_options = optimset('Display', 'final');
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, Erev_lch, iL_corr, S_PCB, on_PCB_cathode, ...
        on_PCB_anode, temp);
    [E_corr, fzval] = fsolve(cor_solver, j0, options);
    i_BV_corr = i_BV(E_corr-Erev_lch, i0, iL_corr, alphas, z, temp);
    i_corr_an = on_PCB_anode.*subplus(i_BV_corr);
    i_corr_cat = on_PCB_cathode.*(-subplus(-i_BV_corr));
    i_corr = i_corr_an + i_corr_cat;
    S_corr = [S_PCB(2:5) sum(S_PCB(2:9)) S_PCB(6:9) sum(S_PCB(2:9)) sum(S_PCB(2:9))];
    I_corr = S_corr.*i_corr;
    
    %%%Solving for mass balance differentials%%%%
    %Calculate concentration/mass balances
    dt = zeros(size(Cm));
    %Electrowinning cell concentration balances
    dt(1) = ((Cm(13)-Cm(1))*Q + I_cell(1)/F/z(1))/vol_cell; %Cu2+
    dt(2) = ((Cm(14)-Cm(2))*Q + I_cell(2)/F/z(2))/vol_cell; %Sn2+
    dt(3) = ((Cm(15)-Cm(3))*Q + I_cell(3)/F/z(3))/vol_cell; %Al3+
    dt(4) = ((Cm(16)-Cm(4))*Q + I_cell(4)/F/z(4))/vol_cell; %Pb2+
    dt(5) = ((Cm(17)-Cm(5))*Q + I_cell(6)/F/z(6)-I_cell(5)/F/z(5))/vol_cell; %Fe2+
    dt(6) = ((Cm(18)-Cm(6))*Q + I_cell(5)/F/z(6))/vol_cell; %Fe3+
    dt(7) = ((Cm(19)-Cm(7))*Q + I_cell(7)/F/z(7))/vol_cell; %Ag+
    dt(8) = ((Cm(20)-Cm(8))*Q + I_cell(8)/F/z(8))/vol_cell; %Au+
    dt(9) = ((Cm(21)-Cm(9))*Q + I_cell(9)/F/z(9))/vol_cell; %Pd2+
    dt(10) = ((Cm(22)-Cm(10))*Q + I_cell(10)/F/z(10) + I_cell(11)/F/z(11))/vol_cell; %H+
    dt(11) = (Cm(23)-Cm(11))*Q/vol_cell; %Cl-
    dt(12) = ((Cm(24)-Cm(12))*Q - 2*I_cell(7)/F/z(7) - 2*I_cell(8)/F/z(8) - 4*I_cell(9)/F/z(9))/vol_cell; %S2O3-
    
    %Leaching vessel concentration balances
    dt(13) = ((Cm(1)-Cm(13))*Q + I_corr(1)/F/z(1))/vol_lch; %Cu2+
    dt(14) = ((Cm(2)-Cm(14))*Q + I_corr(2)/F/z(2))/vol_lch; %Sn2+
    dt(15) = ((Cm(3)-Cm(15))*Q + I_corr(3)/F/z(3))/vol_lch; %Al3+
    dt(16) = ((Cm(4)-Cm(16))*Q + I_corr(4)/F/z(4))/vol_lch; %Pb2+
    dt(17) = ((Cm(5)-Cm(17))*Q + I_corr(6)/F/z(6)-I_corr(5)/F/z(5))/vol_lch; %Fe2+
    dt(18) = ((Cm(6)-Cm(18))*Q + I_corr(5)/F/z(6))/vol_lch; %Fe3+
    dt(19) = ((Cm(7)-Cm(19))*Q + I_corr(7)/F/z(7))/vol_lch; %Ag+
    dt(20) = ((Cm(8)-Cm(20))*Q + I_corr(8)/F/z(8))/vol_lch; %Au+
    dt(21) = ((Cm(9)-Cm(21))*Q + I_corr(9)/F/z(9))/vol_lch; %Pd2+
    dt(22) = ((Cm(10)-Cm(22))*Q + I_corr(10)/F/z(10) + I_corr(11)/F/z(11))/vol_lch; %H+
    dt(23) = (Cm(11)-Cm(23))*Q/vol_lch; %Cl-
    dt(24) = ((Cm(12)-Cm(24))*Q - 2*I_corr(7)/F/z(7) - 2*I_corr(8)/F/z(8) - 4*I_corr(9)/F/z(9))/vol_cell; %S2O3-
    
    %PCB metal mass balances in kg units
    dt(25) = 0; %Inert material 
    dt(26:33) = -mw(2:9).*[I_corr(1:4) I_corr(6:9)]/F./[z(1:4) z(6:9)]/1000; %Solid metals Cu to Pd
    %display current time step
    t
end

function func = cor(Ecorr, Erev, iL, S_PCB, on_PCB_cat, on_PCB_an, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alphas z
    i_BV_corr = i_BV(Ecorr-Erev, i0, iL, alphas, z, temp);
    i_corr_an = on_PCB_an.*subplus(i_BV_corr);
    i_corr_cat = on_PCB_cat.*(-subplus(-i_BV_corr));
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
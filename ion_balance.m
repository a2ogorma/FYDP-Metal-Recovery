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
    
    %universal constants
    global R F;
    R = 8.314; %J/(mol K)
    F = 96485.3329; %C/mol
    
    e = realmin;
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
    %}
    
    %Electrochemical Constants
    % # of electrons in rxn
    global z;
    z_Cu = 2;
    z_Sn = 2;
    z_Al = 3;
    z_Pb = 3;
    z_Fe1 = 1;
    z_Fe2 = 2;
    z_Ag = 1;
    z_Au = 1;
    z_Pd = 2;
    z_H = 0.5;
    z = [z_Cu z_Sn z_Al z_Pb z_Fe1 z_Fe2 z_Ag z_Au z_Pd z_H];
    
    % exchange current densities
    global i0
    i0_Cu = 5E-5; %A/cm^2
    i0_Sn = 1E-6; %A/cm^2
    i0_Al = 1E-6; %A/cm^2
    i0_Pb = 1E-6; %A/cm^2
    i0_Fe1 = 5E-6; %A/cm^2
    i0_Fe2 = 5E-6; %A/cm^2
    i0_Ag = 1E-7; %A/cm^2
    i0_Au = 1E-7; %A/cm^2
    i0_Pd = 1E-8; %A/cm^2
    i0_H = 1E-7; %A/cm^2
    i0 = [i0_Cu i0_Sn i0_Al i0_Pb i0_Fe1 i0_Fe2 i0_Ag i0_Au i0_Pd i0_H];
    
    % charge transfer coefficients
    global alpha
    alpha_Cu = 0.5;
    alpha_Sn = 0.5;
    alpha_Al = 0.5;
    alpha_Pb = 0.5;
    alpha_Fe1 = 0.5;
    alpha_Fe2 = 0.5;
    alpha_Ag = 0.5;
    alpha_Au = 0.5;
    alpha_Pd = 0.5;
    alpha_H = 0.5;
    alpha = [alpha_Cu alpha_Sn alpha_Al alpha_Pb alpha_Fe1 alpha_Fe2...
        alpha_Ag alpha_Au alpha_Pd alpha_H];
    
    % Standard half reaction potentials vs. SHE @ 298 K, 1 atm, 1 M https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    global Eo
    Eo_Cu = 0.337; %V
    Eo_Sn = -0.13; %V
    Eo_Al = -1.66; %V
    Eo_Pb = -0.126; %V
    Eo_Fe1 = 0.77; %V
    Eo_Fe2 = -0.44; %V
    Eo_Ag = 0.7896; %V
    Eo_Au = 1.83; %V
    Eo_Pd = 0.915; %V
    Eo_H = 0; %V
    Eo = [Eo_Cu Eo_Sn Eo_Al Eo_Pb Eo_Fe1 Eo_Fe2 Eo_Ag Eo_Au Eo_Pd Eo_H];
    
    %activity coefficients of ions in solution
    global gamma
    gamma_Cu2 = 1;
    gamma_Fe2 = 1;
    gamma_Fe3 = 1;
    gamma_Sn2 = 1;
    gamma_Al3 = 1;
    gamma_Pb2 = 1;
    gamma_Ag = 1;
    gamma_Au = 1;
    gamma_Pd2 = 1;
    gamma_H = 1;
    gamma_Cl = 1;
    gamma = [gamma_Cu2 gamma_Sn2 gamma_Al3 gamma_Pb2 gamma_Fe2 gamma_Fe3...
        gamma_Ag gamma_Au gamma_Pd2 gamma_H];
    aH2 = 0.0001; %atm
    
    %Lamda infinity values NEED source for iron, tin, nickel rest are from ChE 331 notes
    lamda_Cu2 = 107.2; %S m^2/mol
    lamda_H = 349.81; %S m^2/mol
    lamda_Cl = 76.35; %S m^2/mol
    lamda_Fe2 = 100; %S m^2/mol
    lamda_Fe3 = 100; %S m^2/mol
    lamda_Sn2 = 100; %S m^2/mol
    lamda_Al3 = 100; %S m^2/mol
    lamda_Pb2 = 100; %S m^2/mol
    lamda_Ag = 100; %S m^2/mol
    lamda_Au = 100; %S m^2/mol
    lamda_Pd = 100; %S m^2/mol
    
    %Density of inert material SUBJECT TO CHANGE (Based off range of epoxy)
    rho_inert = 1200; %kg/m^3
    %metal densities in kg/m^3 
    rho_Cu = 8960;
    rho_Sn = 5750;
    rho_Al = 2700;
    rho_Pb = 11350;
    rho_Fe = 7847;
    %rho_Zn = 7133;
    %rho_Ca = 1550;
    %rho_Ni = 8908;
    rho_Ag = 10490;
    rho_Au = 19300;
    rho_Pd = 11900;
    rho = [rho_inert rho_Cu rho_Sn rho_Al rho_Pb rho_Fe rho_Ag rho_Au rho_Pd];
    %molecular weights
    mw_Cu = 64.546;
    mw_Sn = 118.71;
    mw_Al = 26.98;
    mw_Pb = 207.2;
    mw_Fe = 55.845;
    %mw_Zn = 65.38;
    %mw_Ca = 40.078;
    %mw_Ni = 58.6934;
    mw_Ag = 107.8682;
    mw_Au = 196.966;
    mw_Pd = 106.42;
    %index 1 empty to mark place of inert material in other arrays
    mw = [0 mw_Cu mw_Sn mw_Al mw_Pb mw_Fe mw_Ag mw_Au mw_Pd];
    
    %calculate total mass and wt fractions based on partial masses at
    %timestep
    m_PCB = Cm(23:31).';
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
    S_PCB = S_PCB*100^2; %convert m^2 to cm^2
    
    %Nernst potentials - electrowinning cell
    Erev_Cu_cell = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma(1)*max(Cm(1),e)));
    Erev_Sn_cell = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma(2)*max(Cm(2),e)));
    Erev_Al_cell = Eo_Al - R*temp/(z_Al*F)*log(1/(gamma(3)*max(Cm(3),e)));
    Erev_Pb_cell = Eo_Pb - R*temp/(z_Pb*F)*log(1/(gamma(4)*max(Cm(4),e)));
    Erev_Fe1_cell = Eo_Fe1 - R*temp/(z_Fe1*F)*log(gamma(5)*max(Cm(5),e)/(gamma(6)*max(Cm(6),e)));
    Erev_Fe2_cell = Eo_Fe2 - R*temp/(z_Fe2*F)*log(1/(gamma_Fe2*max(Cm(5),e)));
    Erev_Ag_cell = Eo_Ag - R*temp/(z_Ag*F)*log(1/(gamma_Ag*max(Cm(7),e)));
    Erev_Au_cell = Eo_Au - R*temp/(z_Au*F)*log(1/(gamma_Au*max(Cm(8),e)));
    Erev_Pd_cell = Eo_Pd - R*temp/(z_Pd*F)*log(1/(gamma_Pd2*max(Cm(9),e)));
    Erev_H_cell = Eo_H - R*temp/(z_H*F)*log(aH2/(gamma_H*max(Cm(10),e))^2);
    Erev_cell = [Erev_Cu_cell Erev_Sn_cell Erev_Al_cell Erev_Pb_cell Erev_Fe1_cell...
        Erev_Fe2_cell Erev_Ag_cell Erev_Au_cell Erev_Pd_cell Erev_H_cell];
    
    %Nernst potentials - leaching vessel
    Erev_Cu_lch = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(Cm(12),e)));
    Erev_Sn_lch = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma_Sn2*max(Cm(13),e)));
    Erev_Al_lch = Eo_Al - R*temp/(z_Al*F)*log(1/(gamma_Al3*max(Cm(14),e)));
    Erev_Pb_lch = Eo_Pb - R*temp/(z_Pb*F)*log(1/(gamma_Pb2*max(Cm(15),e)));
    Erev_Fe1_lch = Eo_Fe1 - R*temp/(z_Fe1*F)*log(gamma_Fe2*max(Cm(16),e)/gamma_Fe3*max(Cm(17),e));
    Erev_Fe2_lch = Eo_Fe2 - R*temp/(z_Fe2*F)*log(1/(gamma_Fe2*max(Cm(16),e)));
    Erev_Ag_lch = Eo_Ag - R*temp/(z_Ag*F)*log(1/(gamma_Ag*max(Cm(18),e)));
    Erev_Au_lch = Eo_Au - R*temp/(z_Au*F)*log(1/(gamma_Au*max(Cm(19),e)));
    Erev_Pd_lch = Eo_Pd - R*temp/(z_Pd*F)*log(1/(gamma_Pd2*max(Cm(20),e)));
    Erev_H_lch = Eo_H - R*temp/(z_H*F)*log(aH2/(gamma_H*max(Cm(21),e))^2);
    Erev_lch = [Erev_Cu_lch Erev_Sn_lch Erev_Al_lch Erev_Pb_lch Erev_Fe1_lch...
        Erev_Fe2_lch Erev_Ag_lch Erev_Au_lch Erev_Pd_lch Erev_H_lch];
    
    %resistance calculation for IR drop in cell
    kappa = 1000*(Cm(1)*lamda_Cu2+Cm(2)*lamda_Sn2+Cm(3)*lamda_Al3+...
        Cm(4)*lamda_Pb2+Cm(5)*lamda_Fe2+Cm(6)*lamda_Fe3+Cm(7)*lamda_Ag+...
        Cm(8)*lamda_Au+Cm(9)*lamda_Pd+Cm(10)*lamda_H+Cm(11)*lamda_Cl);
    r_sol = l/A_cell/kappa*100;
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol, r_hardware,...
        Erev_cell, S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [0.2, 0.1, -0.1];
    options = optimset('Display','on');
    x = fsolve(solver, x0, options);
    I = x(1);
    E_an = x(2);
    E_cat = x(3);
    eta_cat = E_cat - Erev_cell;
    eta_an = E_an - Erev_cell;
    i_cat = -subplus(-i_BV(eta_cat, i0, alpha, z, temp));
    I_cat = i_cat*S_cat;
    i_an = subplus(i_BV(eta_an(5), i0(5), alpha(5), z(5), temp));
    I_an = i_an*S_an;
    I_cell = I_cat; %overall current for rxn i in cell
    I_cell(5) = I_cell(5)+I_an; 
    Erev_lch
    %solve extraction lch corrosion rate
    j0 = Erev_lch(5)*0.90; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, Erev_lch, S_PCB, temp);
    E_corr = fzero(cor_solver, j0);
    i_corr = i_BV(E_corr-Erev_lch, i0, alpha, z, temp);
    S_corr = [S_PCB(2:5) sum(S_PCB) S_PCB(6:9) 0];
    I_corr = S_corr.*i_corr;
    %Calculate concentration/mass balances
    dt = zeros(size(Cm));
    %Electrowinning cell concentration balances
    dt(1) = ((Cm(12)-Cm(1))*Q + I_cell(1)/F/z(1))/vol_cell; %Cu2+
    dt(2) = ((Cm(13)-Cm(2))*Q + I_cell(2)/F/z(2))/vol_cell; %Sn2+
    dt(3) = ((Cm(14)-Cm(3))*Q + I_cell(3)/F/z(3))/vol_cell; %Al3+
    dt(4) = ((Cm(15)-Cm(4))*Q + I_cell(4)/F/z(4))/vol_cell; %Pb2+
    dt(5) = ((Cm(16)-Cm(5))*Q + I_cell(6)/F/z(6)-I_cell(5)/F/z(5))/vol_cell; %Fe2+
    dt(6) = ((Cm(17)-Cm(6))*Q + I_cell(5)/F/z(6))/vol_cell; %Fe3+
    dt(7) = ((Cm(18)-Cm(7))*Q + I_cell(7)/F/z(7))/vol_cell; %Ag+
    dt(8) = ((Cm(19)-Cm(8))*Q + I_cell(8)/F/z(8))/vol_cell; %Au+
    dt(9) = ((Cm(20)-Cm(9))*Q + I_cell(9)/F/z(9))/vol_cell; %Pd2+
    dt(10) = ((Cm(21)-Cm(10))*Q + I_cell(10)/F/z(10))/vol_lch; %H+
    dt(11) = (Cm(22)-Cm(11))*Q/vol_cell; %Cl-
    
    %Leaching vessel concentration balances
    dt(12) = ((Cm(1)-Cm(12))*Q + I_corr(1)/F/z(1))/vol_lch; %Cu2+
    dt(13) = ((Cm(2)-Cm(13))*Q + I_corr(2)/F/z(2))/vol_lch; %Sn2+
    dt(14) = ((Cm(3)-Cm(14))*Q + I_corr(3)/F/z(3))/vol_lch; %Al3+
    dt(15) = ((Cm(4)-Cm(15))*Q + I_corr(4)/F/z(4))/vol_lch; %Pb2+
    dt(16) = ((Cm(5)-Cm(16))*Q + I_corr(6)/F/z(6)-I_corr(5)/F/z(5))/vol_lch; %Fe2+
    dt(17) = ((Cm(6)-Cm(17))*Q + I_corr(5)/F/z(6))/vol_lch; %Fe3+
    dt(18) = ((Cm(7)-Cm(18))*Q + I_corr(7)/F/z(7))/vol_lch; %Ag+
    dt(19) = ((Cm(8)-Cm(19))*Q + I_corr(8)/F/z(8))/vol_lch; %Au+
    dt(20) = ((Cm(9)-Cm(20))*Q + I_corr(9)/F/z(9))/vol_lch; %Pd2+
    dt(21) = ((Cm(10)-Cm(21))*Q +I_corr(10)/F/z(10))/vol_lch; %H+
    dt(22) = (Cm(11)-Cm(22))*Q/vol_lch; %Cl-
    
    %PCB metal mass balances
    dt(23) = 0; %Inert material 
    dt(24:31) = -mw(2:9).*[I_corr(1:4) I_corr(6:9)]/F./[z(1:4) z(6:9)]; %Solid metals Cu to Pd
    %display current time step
    t
end

function Y = cell_solver(I, E_an, E_cat, V_app, r_sol, r_hardware, Erev, S_an, S_cat, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alpha z
    eta_cat = E_cat - [Erev(1:4) Erev(6:10)];
    eta_an = E_an - Erev(5);
    i_cat = -subplus(-i_BV(eta_cat, [i0(1:4) i0(6:10)], [alpha(1:4) alpha(6:10)], [z(1:4) z(6:10)], temp));
    I_cat = i_cat*S_cat;
    if sum(I_cat) == -Inf
        error("Cathodic current infinite");
    end
    i_an = i_BV(eta_an, i0(5), alpha(5), z(5), temp);
    I_an = i_an*S_an;
    if sum(I_an) == Inf
        error("Anodic current infinite");
    end
    Y(1) = E_an - E_cat + I*(r_sol+r_hardware) - V_app;
    Y(2) = I - I_an;
    Y(3) = I_an + sum(I_cat);
end
%{
function Y = cell_solver(I, E_an, E_cat, V_app, r_sol, r_hardware, Erev, S_an, S_cat, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alpha z
    eta_cat = E_cat - Erev;
    eta_an = E_an - Erev;
    i_cat = -subplus(-i_BV(eta_cat, i0, alpha, z, temp));
    I_cat = i_cat*S_cat;
    i_an = subplus(i_BV(eta_an, i0, alpha, z, temp));
    I_an = i_an*S_an;
    Y(1) = E_an - E_cat + I*(r_sol+r_hardware) - V_app;
    Y(2) = I - sum(I_an);
    Y(3) = sum(I_an) + sum(I_cat);
end
%}

function func = cor(Ecorr, Erev, S_PCB, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alpha z
    i_corr = i_BV(Ecorr-Erev, i0, alpha, z, temp);
    %arrange surface areas in proper order
    S_corr = [S_PCB(2:5) sum(S_PCB) S_PCB(6:9) 0];
    if abs(sum(S_corr.*i_corr)) == Inf
        error("Corrosion currents infinite");
    end
    %imag(sum(S_corr.*i_corr))
    func = sum(S_corr.*i_corr);
end
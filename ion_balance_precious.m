function dC_dt = ion_balance(t, C, temp, pres, vol_cell, vol_bed, Q, S_an, S_cat, V_app, r_particles, l, A_cell)
    %t = time span (s)
    %C = ionic concentration array [C_Cu2+_cell C_Fe2+_cell C_Fe3+_cell C_H+_cell C_Cl-_cell 
    %C_Cu2+_bed C_Fe2+_bed C_Fe3+_bed C_H+_bed C_Cl-_bed] (mol/L)
    
    %universal constants
    global R F;
    R = 8.314; %J/(mol K)
    F = 96485.3329; %C/mol
    %{
        Reactions
        Anodic: Fe2+ --> Fe3+ + e-
        Cathodic: Cu2+ + 2e- --> Cu(s)
    %}
    %Electrochemical Constants
    % # of electrons in rxn
    global z_Au z_Fe z_Ag z_Pd z_An z_H
    z_Au = 1;
    z_Fe = 1;
    z_Ag = 1;
    z_Pd = 2;
    z_An = 4;
    z_H = 1;
    % exchange current densities
    global i0_Au i0_Fe i0_Ag i0_Pd i0_An i0_H
    i0_An = 1E-12; %A/m2
    i0_Fe = 5E-6;%from that reference
    i0_Pd = 7.3e-9; %A/m2, from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1007/s10800-007-9434-x - verify this link applies
    i0_Au = 2E-8; %A/m2, value of dissolution including Na2S in dissolution reaction from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1134/S1023193506040021, converting A/cm2 to A/m2 
    i0_Ag = 4e-9; %using a mix of a few from that reference
    i0_H = 1E-10; %need to identify suitable electrode to help limit hydrogen evolution. From volcano plot in class, likely copper, gold or silver electrodes

    % charge transfer coefficients
    global alpha_Au alpha_Fe alpha_Ag alpha_Pd alpha_An alpha_H
    alpha_Au = 0.5;
    alpha_Fe = 0.5;
    alpha_Ag = 0.5;
    alpha_Pd = 0.5;
    alpha_An = 0.5;
    alpha_H = 0.5;
    % Standard half reaction potentials vs. SHE @ 298 K, 1 atm, 1 M 
    global Eo_Au Eo_Fe Eo_Ag Eo_Pd Eo_An Eo_H
    Eo_Au = 0.153; %V
    Eo_Fe = 0.77; %V
    Eo_Ag = 0.060113; %V
    Eo_Pd = 0.0862; %V
    Eo_An = 1.23; %V
    Eo_H = 0; %V
    %activity coefficients of ions in solution%assume all ideal unless proven otherwise
    global gamma_Au gamma_Fe2 gamma_Fe3 gamma_Ag gamma_Pd gamma_OH gamma_H gamma_S2O3
    gamma_Au = 1;
    gamma_Fe2 = 1;
    gamma_Fe3 = 1;
    gamma_Ag = 1;
    gamma_Pd = 1;
    gamma_OH = 1;
    gamma_H = 1;
    gamma_S2O3 = 1;
    aH2 = 0.0001; %atm
    aO2 = 0.21;%atm

    %Lamda infinity values NEED source for iron, tin, nickel rest are from ChE 331 notes
    lamda_H = 349.81; %S m^2/mol
    lamda_Fe2 = 100; %S m^2/mol
    lamda_Fe3 = 100; %S m^2/mol
    lamda_S2O3 = 100; %S m^2/mol
    lamda_Ag = 61.90e-4;
	lamda_Au = 4.1e7;
	lamda_Pd = 0.1;
    lamda_OH = 0;%easy update
    
    %Surface area calculation for corrosion
    SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
    packing_density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
    S_corr = vol_bed*packing_density*SSA; %m2 Total surface area of spheres in bed.
    
    %Nernst potentials
    Erev_Au_cell = Eo_Au - R*temp/(z_Au*F)*log((gamma_S2O3*max(C(6),eps))^2/(gamma_Au*max(C(1),eps)));
    Erev_Fe_cell = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(2),eps)/gamma_Fe3*max(C(3),eps));
    Erev_Ag_cell = Eo_Ag - R*temp/(z_Ag*F)*log((gamma_S2O3*max(C(6),eps))^2/(gamma_Ag*max(C(4),eps)));
    Erev_Pd_cell = Eo_Pd - R*temp/(z_Pd*F)*log((gamma_S2O3*max(C(6),eps))^4/(gamma_Pd*max(C(5),eps)));
    Erev_An_cell = Eo_An - R*temp/(z_An*F)*log((aO2*(gamma_H*max(C(7),eps))^4));
    Erev_H_cell = Eo_H - R*temp/(z_H*F)*log(aH2^0.5/(gamma_H*max(C(7),eps)));
    
    Erev_Au_bed = Eo_Au - R*temp/(z_Au*F)*log((gamma_S2O3*max(C(13),eps))^2/(gamma_Au*max(C(8),eps)));
    Erev_Fe_bed = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(9),eps)/gamma_Fe3*max(C(10),eps));
    Erev_Ag_bed = Eo_Ag - R*temp/(z_Ag*F)*log((gamma_S2O3*max(C(13),eps))^2/(gamma_Ag*max(C(11),eps)));
    Erev_Pd_bed = Eo_Pd - R*temp/(z_Pd*F)*log((gamma_S2O3*max(C(13),eps))^4/(gamma_Pd*max(C(12),eps)));
    Erev_An_bed = Eo_An - R*temp/(z_An*F)*log((aO2*(gamma_H*max(C(14),eps))^4));
    
    %resistance calculation for IR drops
    kappa = 1000*(C(1)*lamda_Au + C(2)*lamda_Fe2 + C(3)*lamda_Fe3 + C(6)*lamda_S2O3 + C(5)*lamda_Pd + C(7)*lamda_H + C(4)*lamda_Ag);
    r_sol = l/A_cell/kappa*100;
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol, r_hardware, [Erev_Au_cell Erev_Fe_cell Erev_Ag_cell Erev_Pd_cell Erev_An_cell Erev_H_cell], S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    options = optimset('Display','off');
    x = fsolve(solver, x0, options);
    I_an = x(1);
    E_an = x(2);
    E_cat = x(3);
    I_Fe = i_BV(E_an - Erev_Fe_cell, i0_Fe, alpha_Fe, z_Fe, temp);
    I_An = i_BV(E_an - Erev_An_cell, i0_An, alpha_An, z_An, temp);
    I_Au = i_BV(E_cat - Erev_Au_cell, i0_Au, alpha_Au, z_Au, temp);
    I_Ag = i_BV(E_cat - Erev_Ag_cell, i0_Ag, alpha_Ag, z_Ag, temp);
    I_Pd = i_BV(E_cat - Erev_Pd_cell, i0_Pd, alpha_Pd, z_Pd, temp);
    I_H = i_BV(E_an - Erev_H_cell, i0_H, alpha_H, z_H, temp);
    
    %solve extraction bed corrosion rate
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, [Erev_Au_bed, Erev_Fe_bed, Erev_Ag_bed, Erev_Pd_bed Erev_An_bed], temp);
    E_corr = fzero(cor_solver, j0);
    I_corr_Fe = -S_corr*i_BV(E_corr-Erev_Fe_bed, i0_Fe, alpha_Fe, z_Fe, temp);
    I_corr_Au = S_corr*i_BV(E_corr-Erev_Au_bed, i0_Au, alpha_Au, z_Au, temp);
    I_corr_Ag = S_corr*i_BV(E_corr-Erev_Ag_bed, i0_Ag, alpha_Ag, z_Ag, temp);
    I_corr_Pd = S_corr*i_BV(E_corr-Erev_Pd_bed, i0_Pd, alpha_Pd, z_Pd, temp);
    I_corr_An = -S_corr*i_BV(E_corr-Erev_An_bed, i0_An, alpha_An, z_An, temp);
    
    %Calculate concentration balances
    dC_dt = zeros(size(C));
    dC_dt(1) = ((C(8)-C(1))*Q + I_Au/F/z_Au)/vol_cell;
    dC_dt(2) = ((C(9)-C(2))*Q - I_Fe/F/z_Fe)/vol_cell;
    dC_dt(3) = ((C(10)-C(3))*Q + I_Fe/F/z_Fe)/vol_cell;
    dC_dt(4) = ((C(11)-C(4))*Q + I_Ag/F/z_Ag)/vol_cell;
    dC_dt(5) = ((C(12)-C(5))*Q + I_Pd/F/z_Pd)/vol_cell;
    dC_dt(6) = ((C(13)-C(6))*Q - 2*I_Au/F/z_Au - 2*I_Ag/F/z_Ag - 4*I_Pd/F/z_Pd)/vol_cell;
    dC_dt(7) = ((C(14)-C(7))*Q + 4*I_An/F/z_An + I_H/F/z_H)/vol_cell;
    dC_dt(8) = ((C(1)-C(8))*Q + I_corr_Au/F/z_Au)/vol_bed;
    dC_dt(9) = ((C(2)-C(9))*Q + I_corr_Fe/F/z_Fe)/vol_bed;
    dC_dt(10) = ((C(3)-C(10))*Q - I_corr_Fe/F/z_Fe)/vol_bed;
    dC_dt(11) = ((C(4)-C(11))*Q + I_corr_Ag/F/z_Ag)/vol_bed;
    dC_dt(12) = ((C(5)-C(12))*Q + I_corr_Pd/F/z_Pd)/vol_bed;
    dC_dt(13) = ((C(6)-C(13))*Q - 2*I_corr_Au/F/z_Au - 2*I_corr_Ag/F/z_Ag - 4*I_corr_Pd/F/z_Pd)/vol_bed;
    dC_dt(14) = ((C(7)-C(14))*Q + 4*I_An/F/z_An)/vol_bed;
  
    %display current time step
    t
end

function Y = cell_solver(I_an, E_an, E_cat, V_app, r_sol, r_hardware, Erev, S_an, S_cat, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials - 1 = Au, 2 = Fe, 3 = Ag, 4 = Pd, 5 =
    % O2 evolution, 6 = H2 evolution
    global i0_Fe alpha_Fe z_Fe i0_Au alpha_Au z_Au i0_Ag alpha_Ag z_Ag i0_Pd alpha_Pd z_Pd i0_An alpha_An z_An i0_H alpha_H z_H
    eta_Au = E_cat - Erev(1);
    eta_Fe = E_an - Erev(2);
    eta_Ag = E_cat - Erev(3);
    eta_Pd = E_cat - Erev(4);
    eta_An = E_an - Erev(5);
    eta_H = E_cat - Erev(6);
    i_cat = i_BV(eta_Au, i0_Au, alpha_Au, z_Au, temp) + i_BV(eta_Ag, i0_Ag, alpha_Ag, z_Ag, temp) + i_BV(eta_Pd, i0_Pd, alpha_Pd, z_Pd, temp) + i_BV(eta_H, i0_H, alpha_H, z_H, temp);
    Y(1) = E_an - E_cat + I_an*(r_sol+r_hardware) - V_app;
    Y(2) = I_an - (i_BV(eta_Fe, i0_Fe, alpha_Fe, z_Fe, temp)+i_BV(eta_An, i0_An, alpha_An, z_An, temp))*S_an;
    Y(3) = I_an + i_cat*S_cat;
end

function [func] = cor(Ecorr, Erev, temp)
  %1 = Cu, 2= Fe, 3= Sn, 4 = Ni
  global i0_Fe alpha_Fe z_Fe i0_Au alpha_Au z_Au i0_Ag alpha_Ag z_Ag i0_Pd alpha_Pd z_Pd i0_An alpha_An z_An
  i_Au = i_BV(Ecorr-Erev(1), i0_Au, alpha_Au, z_Au, temp);
  i_Fe = i_BV(Ecorr-Erev(2), i0_Fe, alpha_Fe, z_Fe, temp);
  i_Ag = i_BV(Ecorr-Erev(3), i0_Ag, alpha_Ag, z_Ag, temp);
  i_Pd = i_BV(Ecorr-Erev(4), i0_Pd, alpha_Pd, z_Pd, temp);
  i_An = i_BV(Ecorr-Erev(5), i0_An, alpha_An, z_An, temp);
  func = i_Au + i_Fe + i_Ag + i_Pd + i_An;
  end
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
    global z_Cu z_Fe z_Sn z_Ni 
    z_Cu = 2;
    z_Fe = 1;
    z_Sn = 2;
    z_Ni = 2;
    % exchange current densities
    global i0_Cu i0_Fe i0_Sn i0_Ni
    i0_Cu = 5E-5; %A/cm^2
    i0_Fe = 5E-6; %A/cm^2
    i0_Sn = 1E-6; %A/cm^2
    i0_Ni = 1E-8; %A/cm^2
    % charge transfer coefficients
    global alpha_Cu alpha_Fe alpha_Sn alpha_Ni
    alpha_Cu = 0.5;
    alpha_Fe = 0.5;
    alpha_Sn = 0.5;
    alpha_Ni = 0.5;
    % Standard half reaction potentials vs. SHE @ 298 K, 1 atm, 1 M https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    global Eo_Cu Eo_Fe Eo_Sn Eo_Ni
    Eo_Cu = 0.337; %V
    Eo_Fe = 0.77; %V
    Eo_Sn = -0.13; %V
    Eo_Ni = -0.25; %V
    %activity coefficients of ions in solution
    global gamma_Cu2 gamma_Fe2 gamma_Fe3 gamma_Sn2 gamma_Ni2
    gamma_Cu2 = 1;
    gamma_Fe2 = 1;
    gamma_Fe3 = 1;
    gamma_Sn2 = 1;
    gamma_Ni2 = 1;
    %Lamda infinity values NEED source for iron, tin, nickel rest are from ChE 331 notes
    lamda_Cu2 = 107.2; %S m^2/mol
    lamda_H = 349.81; %S m^2/mol
    lamda_Cl = 76.35; %S m^2/mol
    lamda_Fe2 = 100; %S m^2/mol
    lamda_Fe3 = 100; %S m^2/mol
    lamda_Sn2 = 100; %S m^2/mol
    lamda_Ni2 = 100; %S m^2/mol
    
    %Surface area calculation for corrosion
    SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
    packing_density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
    S_corr = vol_bed*packing_density*SSA; %m2 Total surface area of spheres in bed.
    
    %Nernst potentials
    Erev_Cu_cell = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(1),eps)));
    Erev_Fe_cell = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(2),eps)/gamma_Fe3*max(C(3),eps));
    Erev_Sn_cell = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma_Sn2*max(C(4),eps)));
    Erev_Ni_cell = Eo_Ni - R*temp/(z_Ni*F)*log(1/(gamma_Ni2*max(C(5),eps)));
    Erev_Cu_bed = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(8),eps)));
    Erev_Fe_bed = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(9),eps)/gamma_Fe3*max(C(10),eps));
    Erev_Sn_bed = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma_Sn2*max(C(11),eps)));
    Erev_Ni_bed = Eo_Ni - R*temp/(z_Ni*F)*log(1/(gamma_Ni2*max(C(12),eps)));
    
    %resistance calculation for IR drops
    kappa = 1000*(C(1)*lamda_Cu2 + C(2)*lamda_Fe2 + C(3)*lamda_Fe3 + C(4)*lamda_Sn2 + C(5)*lamda_Ni2 + C(6)*lamda_H + C(7)*lamda_Cl);
    r_sol = l/A_cell/kappa*100;
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol, r_hardware, [Erev_Cu_cell Erev_Fe_cell Erev_Sn_cell Erev_Ni_cell], S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    options = optimset('Display','off');
    x = fsolve(solver, x0, options);
    I_an = x(1);
    E_an = x(2);
    E_cat = x(3);
    I_Fe = I_an;
    I_Cu = i_BV(E_cat - Erev_Cu_cell, i0_Cu, alpha_Cu, z_Cu, temp);
    I_Sn = i_BV(E_cat - Erev_Sn_cell, i0_Sn, alpha_Sn, z_Sn, temp);
    I_Ni = i_BV(E_cat - Erev_Ni_cell, i0_Ni, alpha_Ni, z_Ni, temp);
    
    %solve extraction bed corrosion rate
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, [Erev_Cu_bed, Erev_Fe_bed, Erev_Sn_bed, Erev_Ni_bed], temp);
    E_corr = fzero(cor_solver, j0);
    I_corr = -S_corr*i_BV(E_corr-Erev_Fe_bed, i0_Fe, alpha_Fe, z_Fe, temp);
    I_corr_Cu = S_corr*i_BV(E_corr-Erev_Cu_bed, i0_Cu, alpha_Cu, z_Cu, temp);
    I_corr_Sn = S_corr*i_BV(E_corr-Erev_Sn_bed, i0_Sn, alpha_Sn, z_Sn, temp);
    I_corr_Ni = S_corr*i_BV(E_corr-Erev_Ni_bed, i0_Ni, alpha_Ni, z_Ni, temp);
    
    %Calculate concentration balances
    dC_dt = zeros(size(C));
    dC_dt(1) = ((C(8)-C(1))*Q + I_Cu/F/z_Cu)/vol_cell;
    dC_dt(2) = ((C(9)-C(2))*Q - I_Fe/F/z_Fe)/vol_cell;
    dC_dt(3) = ((C(10)-C(3))*Q + I_Fe/F/z_Fe)/vol_cell;
    dC_dt(4) = ((C(11)-C(4))*Q + I_Sn/F/z_Sn)/vol_cell;
    dC_dt(5) = ((C(12)-C(5))*Q + I_Ni/F/z_Ni)/vol_cell;
    dC_dt(6) = (C(13)-C(6))*Q/vol_cell;
    dC_dt(7) = (C(14)-C(7))*Q/vol_cell;
    dC_dt(8) = ((C(1)-C(8))*Q + I_corr_Cu/F/z_Cu)/vol_bed;
    dC_dt(9) = ((C(2)-C(9))*Q + I_corr/F/z_Fe)/vol_bed;
    dC_dt(10) = ((C(3)-C(10))*Q - I_corr/F/z_Fe)/vol_bed;
    dC_dt(11) = ((C(4)-C(11))*Q + I_corr_Sn/F/z_Sn)/vol_bed;
    dC_dt(12) = ((C(5)-C(12))*Q + I_corr_Ni/F/z_Ni)/vol_bed;
    dC_dt(13) = (C(6)-C(13))*Q/vol_bed;
    dC_dt(14) = (C(7)-C(14))*Q/vol_bed;
  
    %display current time step
    t
end

function Y = cell_solver(I_an, E_an, E_cat, V_app, r_sol, r_hardware, Erev, S_an, S_cat, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials - 1 = Cu, 2 = Fe, 3 = Sn, 4 = Ni
    global i0_Fe alpha_Fe z_Fe i0_Cu alpha_Cu z_Cu i0_Sn alpha_Sn z_Sn i0_Ni alpha_Ni z_Ni
    eta_an = E_an - Erev(1);
    eta_Cu = E_cat - Erev(2);
    eta_Sn = E_cat - Erev(3);
    eta_Ni = E_cat - Erev(4);
    i_cat = i_BV(eta_Cu, i0_Cu, alpha_Cu, z_Cu, temp) + i_BV(eta_Sn, i0_Sn, alpha_Sn, z_Sn, temp) + i_BV(eta_Ni, i0_Ni, alpha_Ni, z_Ni, temp);
    Y(1) = E_an - E_cat + I_an*(r_sol+r_hardware) - V_app;
    Y(2) = I_an - i_BV(eta_an, i0_Fe, alpha_Fe, z_Fe, temp)*S_an;
    Y(3) = I_an + i_cat*S_cat;
end

function [func] = cor(Ecorr, Erev, temp)
  %1 = Cu, 2= Fe, 3= Sn, 4 = Ni
  global i0_Fe alpha_Fe z_Fe i0_Cu alpha_Cu z_Cu i0_Sn alpha_Sn z_Sn i0_Ni alpha_Ni z_Ni
  i_Cu = i_BV(Ecorr-Erev(1), i0_Cu, alpha_Cu, z_Cu, temp);
  i_Fe = i_BV(Ecorr-Erev(2), i0_Fe, alpha_Fe, z_Fe, temp);
  i_Sn = i_BV(Ecorr-Erev(3), i0_Sn, alpha_Sn, z_Sn, temp);
  i_Ni = i_BV(Ecorr-Erev(4), i0_Ni, alpha_Ni, z_Ni, temp);
  func = i_Cu + i_Fe + i_Sn + i_Ni;
  end
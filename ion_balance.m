function dC_dt = ion_balance(t, C, temp, pres, vol_cell, vol_bed, Q, S_an, S_cat, V_app, r_particles, l, A_cell)
    %t = time span (s)
    %C = ionic concentration array [C_Cu2+_cell C_Fe2+_cell C_Fe3+_cell C_H+_cell C_Cl-_cell 
    %C_Cu2+_bed C_Fe2+_bed C_Fe3+_bed C_H+_bed C_Cl-_bed] (mol/L)
    
    %universal constants
    global R F;
    R = 8.314; %J/(mol K)
    F = 96485.3329; %C/mol
    t
    %{
        Reactions
        Anodic: Fe2+ --> Fe3+ + e-
        Cathodic: Cu2+ + 2e- --> Cu(s)
    %}
    %Defining constants as global variables for use in cell solver function
    global z_Cu z_Fe i0_Cu i0_Fe alpha_Cu alpha_Fe Eo_Cu Eo_Fe gamma_Cu2 gamma_Fe2 gamma_Fe3
    % # of electrons in rxn
    z_Cu = 2;
    z_Fe = 1;
    n = z_Cu+z_Fe;
    % exchange current densities
    i0_Cu = 5E-5; %A/cm^2
    i0_Fe = 5E-8; %A/cm^2
    % charge transfer coefficients
    alpha_Cu = 0.5;
    alpha_Fe = 0.5;
    % Standard half reaction potentials vs. SHE @ 298 K, 1 atm, 1 M https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    Eo_Cu = 0.337; %V
    Eo_Fe = 0.77; %V
    %activity coefficients of ions in solution
    gamma_Cu2 = 1;
    gamma_Fe2 = 1;
    gamma_Fe3 = 1;
    %Lamda infinity values NEED source for iron, rest are from 331 notes
    lamda_Cu2 = 107.2; %S m^2/mol
    lamda_H = 349.81; %S m^2/mol
    lamda_Cl = 76.35; %S m^2/mol
    lamda_Fe2 = 100; %S m^2/mol
    lamda_Fe3 = 100; %S m^2/mol
    
    %Surface area calculation for corrosion
    SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
    packing_density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
    S_corr = vol_bed*packing_density*SSA; %m2 Total surface area of spheres in bed.
    
    %Nernst potentials
    Erev_Cu_cell = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(1),eps)));
    Erev_Fe_cell = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(2),eps)/gamma_Fe3*max(C(3),eps));
    Erev_Cu_bed = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(6),eps)));
    Erev_Fe_bed = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(7),eps)/gamma_Fe3*max(C(8),eps));
    
    %resistance calculation for IR drops
    kappa = 1000*(C(1)*lamda_Cu2 + C(2)*lamda_Fe2 + C(3)*lamda_Fe3 + C(4)*lamda_H + C(5)*lamda_Cl);
    r_sol = l/A_cell/kappa*100;
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol, r_hardware, Erev_Fe_cell, Erev_Cu_cell, S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    options = optimset('Display','off');
    x = fsolve(solver, x0, options);
    I_an = x(1);
    E_an = x(2);
    E_cat = x(3);
    I_Cu = - I_an;
    I_Fe = I_an;
    
    %solve extraction bed corrosion rate
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, [Erev_Cu_bed, Erev_Fe_bed], temp);
    E_corr = fzero(cor_solver, j0);
    I_corr = S_corr*i_BV(E_corr-Erev_Cu_bed, i0_Cu, alpha_Cu, z_Cu, temp);
    
    %Calculate concentration balances
    dC_dt = zeros(size(C));
    dC_dt(1) = ((C(6)-C(1))*Q + I_Cu*F/z_Cu)/vol_cell;
    dC_dt(2) = ((C(7)-C(2))*Q - I_Fe*F/z_Fe)/vol_cell;
    dC_dt(3) = ((C(8)-C(3))*Q + I_Fe*F/z_Fe)/vol_cell;
    dC_dt(4) = (C(9)-C(4))*Q/vol_cell;
    dC_dt(5) = (C(10)-C(5))*Q/vol_cell;
    dC_dt(6) = ((C(1)-C(6))*Q + I_corr*F/z_Cu)/vol_bed;
    dC_dt(7) = ((C(2)-C(7))*Q + I_corr*F/z_Fe)/vol_bed;
    dC_dt(8) = ((C(3)-C(8))*Q - I_corr*F/z_Fe)/vol_bed;
    dC_dt(9) = (C(4)-C(9))*Q/vol_bed;
    dC_dt(10) = (C(5)-C(10))*Q/vol_bed;
end

function Y = cell_solver(I_an, E_an, E_cat, V_app, r_sol, r_hardware, Erev_Fe, Erev_Cu, S_an, S_cat, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    global i0_Fe alpha_Fe z_Fe i0_Cu alpha_Cu z_Cu
    eta_an = E_an - Erev_Fe;
    eta_cat = E_cat - Erev_Cu;
    Y(1) = E_an - E_cat + I_an*(r_sol+r_hardware) - V_app;
    Y(2) = I_an - i_BV(eta_an, i0_Fe, alpha_Fe, z_Fe, temp)*S_an;
    Y(3) = I_an + i_BV(eta_cat, i0_Cu, alpha_Cu, z_Cu, temp)*S_cat;
end

function [func] = cor(Ecorr, Erev, temp)
  %1 = Cu, 2 = Fe
  global i0_Cu i0_Fe alpha_Cu alpha_Fe z_Cu z_Fe
  i_Cu = i_BV(Ecorr-Erev(1), i0_Cu, alpha_Cu, z_Cu, temp);
  i_Fe = i_BV(Ecorr-Erev(2), i0_Fe, alpha_Fe, z_Fe, temp);
  func = i_Cu + i_Fe;
  end
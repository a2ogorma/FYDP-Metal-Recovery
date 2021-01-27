%{
	Continuous stirred tank model for base metal extraction/recovery cell
    Using ode45
	To Add: Applied Voltage calc. Non-constant volume?
%} 
clear variables

%universal constants
global R F;
R = 8.314; %J/(mol K)
F = 96485.3329; %C/mol

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

% system parameters
temp = 298; %K
pres = 1; % atm
vol_cell = 250; %L
Q = 5; % L/s (flowrate)
%Electrode areas
S_cat = 500; %cm^2
S_an = 500; %cm^2
%Cross sectional area of cell
A_cell = 500; %cm^2
%Length b/w electrodes
l = 100; %cm
%Applied Voltage (potentiostat)
V_app = 12; %V
%Extraction vessel parameters
vol_bed = 200; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
r_particles = 0.001; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
tfinal = 12; %s

%Surface area calculation for corrosion
SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
packing_density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
S_corr = 0.001*vol_bed*packing_density*SSA; %m2 Total surface area of spheres in bed.

%initial concentrations in mol/L
%Cell Concentrations (recovery)
Ci_Cu2_cell = 0.2;
Ci_Fe2_cell = 0.5;
Ci_Fe3_cell = 0.001;
Ci_Sn2_cell = 0.001;
Ci_Ni2_cell = 0.001;
Ci_H_cell = 0.5;
Ci_Cl_cell = 2*(Ci_Cu2_cell+Ci_Fe2_cell)+Ci_H_cell; %Calculation to ensure 
%electrolyte has net neutral charge

%Bed concentrations (extraction)
Ci_Cu2_bed = eps;%0.2;
Ci_Fe2_bed = 0.5;
Ci_Fe3_bed = 0.001;
Ci_Sn2_bed = 0.001;
Ci_Ni2_bed = 0.001;
Ci_H_bed = 0;
Ci_Cl_bed = 2*(Ci_Cu2_bed+Ci_Fe2_bed)+Ci_H_bed;

%initializing concentrations [C_Cu2+_cell C_Fe2+_cell C_Fe3+_cell C_H+_cell C_Cl-_cell 
%C_Cu2+_bed C_Fe2+_bed C_Fe3+_bed C_H+_bed C_Cl-_bed] (mol/L)
Ci = [Ci_Cu2_cell Ci_Fe2_cell Ci_Fe3_cell Ci_Sn2_cell Ci_Ni2_cell Ci_H_cell Ci_Cl_cell Ci_Cu2_bed Ci_Fe2_bed Ci_Fe3_bed Ci_Sn2_bed Ci_Ni2_bed Ci_H_bed Ci_Cl_bed];

%solve conc profiles
tspan = [0 tfinal];
options = odeset('NonNegative',1:10);
balance_solver = @(t, C) ion_balance(t, C, temp, pres, vol_cell, vol_bed, Q, S_an, S_cat, V_app, r_particles, l, A_cell);
[t, C] = ode15s(balance_solver, tspan, Ci, options);

%backcalculate currents/potentials
Erev_Cu_cell = zeros(size(t));
Erev_Cu_bed = zeros(size(t));
Erev_Fe_cell = zeros(size(t));
Erev_Fe_bed = zeros(size(t));
Erev_Sn_cell = zeros(size(t));
Erev_Sn_bed = zeros(size(t));
Erev_Ni_cell = zeros(size(t));
Erev_Ni_bed = zeros(size(t));
I_an = zeros(size(t));
E_an = zeros(size(t));
E_cat = zeros(size(t));
I_Fe = zeros(size(t));
I_corr = zeros(size(t));
E_corr = zeros(size(t));
r_sol = zeros(size(t));

for j = 1:1:length(t)
    disp(t(j))
    %Nernst Potentials
    Erev_Cu_cell(j) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(j,1),eps)));
    Erev_Fe_cell(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(j,2),eps)/gamma_Fe3*max(C(j,3),eps));
    Erev_Sn_cell(j) = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma_Sn2*max(C(j,4),eps)));
    Erev_Ni_cell(j) = Eo_Ni - R*temp/(z_Ni*F)*log(1/(gamma_Ni2*max(C(j,5),eps)));
    Erev_Cu_bed(j) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(j,8),eps)));
    Erev_Fe_bed(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(j,9),eps)/gamma_Fe3*max(C(j,10),eps));
    Erev_Sn_bed(j) = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma_Sn2*max(C(j,11),eps)));
    Erev_Ni_bed(j) = Eo_Ni - R*temp/(z_Ni*F)*log(1/(gamma_Ni2*max(C(j,12),eps)));
    
    %resistance calculation for IR drops
    kappa = 1000*(C(j,1)*lamda_Cu2 + C(j,2)*lamda_Fe2 + C(j,3)*lamda_Fe3 + C(j,4)*lamda_Sn2 + C(j,5)*lamda_Ni2 + C(j,6)*lamda_H + C(j,7)*lamda_Cl);
    r_sol(j) = l/A_cell/kappa*100; %ohms
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol(j), r_hardware, [Erev_Cu_cell(j) Erev_Fe_cell(j) Erev_Sn_cell(j) Erev_Ni_cell(j)], S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    options = optimset('Display','off');
    x = fsolve(solver, x0, options);
    I_an(j) = x(1);
    E_an(j) = x(2);
    E_cat(j) = x(3);
    I_Fe(j) = I_an(j);
    
    %solve extraction bed corrosion rate
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, [Erev_Cu_bed(j), Erev_Fe_bed(j), Erev_Sn_bed(j), Erev_Ni_bed(j)], temp);
    E_corr(j) = fzero(cor_solver, j0);
    I_corr(j) = -S_corr*i_BV(E_corr(j)-Erev_Fe_bed(j), i0_Fe, alpha_Fe, z_Fe, temp);
end

%Calculate currents, overpotentials
eta_Cu = E_cat - Erev_Cu_cell;
eta_Sn = E_cat - Erev_Sn_cell;
eta_Ni = E_cat - Erev_Ni_cell;
I_Cu = S_cat*i_BV(eta_Cu, i0_Cu, alpha_Cu, z_Cu, temp);
I_Sn = S_cat*i_BV(eta_Sn, i0_Sn, alpha_Sn, z_Sn, temp);
I_Ni = S_cat*i_BV(eta_Cu, i0_Ni, alpha_Ni, z_Ni, temp);

%plots
subplot(3,2,1)
plot(t,C(:,1),t,C(:,8))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Copper')
subplot(3,2,2)
plot(t,C(:,2),t,C(:,9))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(II)')
subplot(3,2,3)
plot(t,C(:,3),t,C(:,10))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(III)')
subplot(3,2,4)
plot(t,C(:,4),t,C(:,11))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Tin')
subplot(3,2,5)
plot(t,C(:,5),t,C(:,12))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Nickel')
subplot(3,2,6)
plot(t,I_corr,t,I_Cu,t,I_Fe)
xlabel('Time (s)')
ylabel('Current (A)')
title('Currents')
legend('Corr','CopperCell','IronCell')

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

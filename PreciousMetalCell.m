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
Ci_Au_cell = eps;
Ci_Fe2_cell = 0.05;
Ci_Fe3_cell = 0.05;
Ci_Ag_cell = eps;
Ci_Pd_cell = eps;
Ci_H_cell = 1e-10;
Ci_S2O3_cell = 0.01;

%Bed concentrations (extraction)
Ci_Au_bed = eps;%0.2;
Ci_Fe2_bed = 0.5;
Ci_Fe3_bed = 0.001;
Ci_Ag_bed = 0.001;
Ci_Pd_bed = 0.001;
Ci_H_bed = 1e-10;
Ci_S2O3_bed = 0.5;

%initializing concentrations [C_Cu2+_cell C_Fe2+_cell C_Fe3+_cell C_H+_cell C_Cl-_cell 
%C_Cu2+_bed C_Fe2+_bed C_Fe3+_bed C_H+_bed C_Cl-_bed] (mol/L)
Ci = [Ci_Au_cell Ci_Fe2_cell Ci_Fe3_cell Ci_Ag_cell Ci_Pd_cell Ci_H_cell Ci_S2O3_cell Ci_Au_bed Ci_Fe2_bed Ci_Fe3_bed Ci_Ag_bed Ci_Pd_bed Ci_H_bed Ci_S2O3_bed];

%solve conc profiles
tspan = [0 tfinal];
options = odeset('NonNegative',1:10);
balance_solver = @(t, C) ion_balance_precious(t, C, temp, pres, vol_cell, vol_bed, Q, S_an, S_cat, V_app, r_particles, l, A_cell);
[t, C] = ode15s(balance_solver, tspan, Ci, options);

%backcalculate currents/potentials
Erev_Au_cell = zeros(size(t));
Erev_Au_bed = zeros(size(t));
Erev_Fe_cell = zeros(size(t));
Erev_Fe_bed = zeros(size(t));
Erev_Ag_cell = zeros(size(t));
Erev_Ag_bed = zeros(size(t));
Erev_Pd_cell = zeros(size(t));
Erev_Pd_bed = zeros(size(t));
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
    Erev_Au_cell(j) = Eo_Au - R*temp/(z_Au*F)*log((gamma_S2O3*max(C(j,6),eps))^2/(gamma_Au*max(C(j,1),eps)));
    Erev_Fe_cell(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(j,2),eps)/gamma_Fe3*max(C(j,3),eps));
    Erev_Ag_cell(j) = Eo_Ag - R*temp/(z_Ag*F)*log((gamma_S2O3*max(C(j,6),eps))^2/(gamma_Ag*max(C(j,4),eps)));
    Erev_Pd_cell(j) = Eo_Pd - R*temp/(z_Pd*F)*log((gamma_S2O3*max(C(j,6),eps))^4/(gamma_Pd*max(C(j,5),eps)));
    Erev_An_cell(j) = Eo_An - R*temp/(z_An*F)*log((aO2*(gamma_H*max(C(j,7),eps))^4));
    Erev_H_cell(j) = Eo_H - R*temp/(z_H*F)*log(aH2^0.5/(gamma_H*max(C(j,7),eps)));
    
    Erev_Au_bed(j) = Eo_Au - R*temp/(z_Au*F)*log((gamma_S2O3*max(C(j,13),eps))^2/(gamma_Au*max(C(j,8),eps)));
    Erev_Fe_bed(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(j,9),eps)/gamma_Fe3*max(C(j,10),eps));
    Erev_Ag_bed(j) = Eo_Ag - R*temp/(z_Ag*F)*log((gamma_S2O3*max(C(j,13),eps))^2/(gamma_Ag*max(C(j,11),eps)));
    Erev_Pd_bed(j) = Eo_Pd - R*temp/(z_Pd*F)*log((gamma_S2O3*max(C(j,13),eps))^4/(gamma_Pd*max(C(j,12),eps)));
    Erev_An_bed(j) = Eo_An - R*temp/(z_An*F)*log((aO2*(gamma_H*max(C(j,14),eps))^4));
    
    %resistance calculation for IR drops
    kappa = 1000*(C(j,1)*lamda_Au + C(j,2)*lamda_Fe2 + C(j,3)*lamda_Fe3 + C(j,6)*lamda_S2O3 + C(j,5)*lamda_Pd + C(j,7)*lamda_H + C(j,4)*lamda_Ag);
    r_sol(j) = l/A_cell/kappa*100; %ohms
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol(j), r_hardware, [Erev_Au_cell(j) Erev_Fe_cell(j) Erev_Ag_cell(j) Erev_Pd_cell(j) Erev_An_cell(j) Erev_H_cell(j)], S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    options = optimset('Display','off');
    x = fsolve(solver, x0, options);
    I_an(j) = x(1);
    E_an(j) = x(2);
    E_cat(j) = x(3);
    I_Fe(j) = i_BV(E_an(j) - Erev_Fe_cell(j), i0_Fe, alpha_Fe, z_Fe, temp);
    I_An(j) = i_BV(E_an(j) - Erev_An_cell(j), i0_An, alpha_An, z_An, temp);
    I_Au(j) = i_BV(E_cat(j) - Erev_Au_cell(j), i0_Au, alpha_Au, z_Au, temp);
    I_Ag(j) = i_BV(E_cat(j) - Erev_Ag_cell(j), i0_Ag, alpha_Ag, z_Ag, temp);
    I_Pd(j) = i_BV(E_cat(j) - Erev_Pd_cell(j), i0_Pd, alpha_Pd, z_Pd, temp);
    I_H(j) = i_BV(E_an(j) - Erev_H_cell(j), i0_H, alpha_H, z_H, temp);
    
    
    %solve extraction bed corrosion rate
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, [Erev_Au_bed(j), Erev_Fe_bed(j), Erev_Ag_bed(j), Erev_Pd_bed(j) Erev_An_bed(j)], temp);
    E_corr(j) = fzero(cor_solver, j0);
    I_corr_Fe(j) = -S_corr*i_BV(E_corr(j)-Erev_Fe_bed(j), i0_Fe, alpha_Fe, z_Fe, temp);
    I_corr_Au(j) = S_corr*i_BV(E_corr(j)-Erev_Au_bed(j), i0_Au, alpha_Au, z_Au, temp);
    I_corr_Ag(j) = S_corr*i_BV(E_corr(j)-Erev_Ag_bed(j), i0_Ag, alpha_Ag, z_Ag, temp);
    I_corr_Pd(j) = S_corr*i_BV(E_corr(j)-Erev_Pd_bed(j), i0_Pd, alpha_Pd, z_Pd, temp);
    I_corr_An(j) = -S_corr*i_BV(E_corr(j)-Erev_An_bed(j), i0_An, alpha_An, z_An, temp);
    
end

%Calculate currents, overpotentials
eta_Cu = E_cat - Erev_Au_cell;
eta_Sn = E_cat - Erev_Ag_cell;
eta_Ni = E_cat - Erev_Pd_cell;

%plots
subplot(3,2,1)
plot(t,C(:,1),t,C(:,8))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Gold')
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
title('Silver')
subplot(3,2,5)
plot(t,C(:,5),t,C(:,12))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Palladium')
subplot(3,2,6)
plot(t,I_corr,t,(I_Au+I_Ag+I_Pd+I_H),t,(I_Fe+I_An))
xlabel('Time (s)')
ylabel('Current (A)')
title('Currents')
legend('Corr','Cathode','Anode')

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
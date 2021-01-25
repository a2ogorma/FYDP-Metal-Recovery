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
global z_Cu z_Fe i0_Cu i0_Fe alpha_Cu alpha_Fe Eo_Cu Eo_Fe gamma_Cu2 gamma_Fe2 gamma_Fe3
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
vol_bed = 0.2; %m3 (Initial) volume of bed holding the particles assuming the bed is completly full.
r_particles = 0.001; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.

%Surface area calculation for corrosion
SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
packing_density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
S_corr = vol_bed*packing_density*SSA; %m2 Total surface area of spheres in bed.

%initial concentrations in mol/L
%Cell Concentrations (recovery)
Ci_Cu2_cell = 0.2;
Ci_Fe2_cell = 0.5;
Ci_Fe3_cell = 0.001;
Ci_H_cell = 0.5;
Ci_Cl_cell = 2*(Ci_Cu2_cell+Ci_Fe2_cell)+Ci_H_cell; %Calculation to ensure 
%electrolyte has net neutral charge

%Bed concentrations (extraction)
Ci_Cu2_bed = 0.2;
Ci_Fe2_bed = 0.5;
Ci_Fe3_bed = 0.001;
Ci_H_bed = 0;
Ci_Cl_bed = 2*(Ci_Cu2_bed+Ci_Fe2_bed)+Ci_H_bed;

%initializing concentrations [C_Cu2+_cell C_Fe2+_cell C_Fe3+_cell C_H+_cell C_Cl-_cell 
%C_Cu2+_bed C_Fe2+_bed C_Fe3+_bed C_H+_bed C_Cl-_bed] (mol/L)
Ci = [Ci_Cu2_cell Ci_Fe2_cell Ci_Fe3_cell Ci_H_cell Ci_Cl_cell Ci_Cu2_bed Ci_Fe2_bed Ci_Fe3_bed Ci_H_bed Ci_Cl_bed];

%solve conc profiles
tspan = [0 50];
balance_solver = @(t, C) ion_balance(t, C, temp, pres, vol_cell, vol_bed, Q, S_an, S_cat, V_app, r_particles, l, A_cell);
[t, C] = ode45(balance_solver, tspan, Ci);

%backcalculate currents/potentials
Erev_Cu_cell = zeros(size(t));
Erev_Cu_bed = zeros(size(t));
Erev_Fe_cell = zeros(size(t));
Erev_Fe_bed = zeros(size(t));
I_an = zeros(size(t));
E_an = zeros(size(t));
E_cat = zeros(size(t));
I_Cu = zeros(size(t));
I_Fe = zeros(size(t));
I_corr = zeros(size(t));
E_corr = zeros(size(t));
r_sol = zeros(size(t));

for j = 1:1:length(t)
    %Nernst Potentials
    Erev_Cu_cell(j) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(j,1),eps)));
    Erev_Fe_cell(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(j,2),eps)/gamma_Fe3*max(C(j,3),eps));
    Erev_Cu_bed(j) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(C(j,6),eps)));
    Erev_Fe_bed(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(C(j,7),eps)/gamma_Fe3*max(C(j,8),eps));
    
    %resistance calculation for IR drops
    kappa = 1000*(C(j,1)*lamda_Cu2 + C(j,2)*lamda_Fe2 + C(j,3)*lamda_Fe3 + C(j,4)*lamda_H + C(j,5)*lamda_Cl);
    r_sol(j) = l/A_cell/kappa*100; %ohms
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol(j), r_hardware, Erev_Fe_cell(j), Erev_Cu_cell(j), S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    options = optimset('Display','off');
    x = fsolve(solver, x0, options);
    I_an(j) = x(1);
    E_an(j) = x(2);
    E_cat(j) = x(3);
    I_Cu(j) = - I_an(j);
    I_Fe(j) = I_an(j);
    
    %solve extraction bed corrosion rate
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, [Erev_Cu_bed(j), Erev_Fe_bed(j)], temp);
    E_corr(j) = fzero(cor_solver, j0);
    I_corr(j) = S_corr*i_BV(E_corr(j)-Erev_Cu_bed(j), i0_Cu, alpha_Cu, z_Cu, temp);
end
eta_an = E_an - Erev_Fe_cell;
eta_cat = E_cat - Erev_Cu_cell;

%plots
subplot(3,2,1)
plot(t,C(:,1),t,C(:,6))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Copper')
subplot(3,2,2)
plot(t,C(:,2),t,C(:,7))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(II)')
subplot(3,2,3)
plot(t,C(:,3),t,C(:,8))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(III)')
subplot(3,2,4)
plot(t,C(:,4),t,C(:,9))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Hydrogen')
subplot(3,2,5)
plot(t,C(:,5),t,C(:,10))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Chloride')
subplot(3,2,6)
plot(t,I_corr)
xlabel('Time (s)')
ylabel('Corrosion Current (A)')
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
  %1 = Cu, 2= Fe
  global i0_Cu i0_Fe alpha_Cu alpha_Fe z_Cu z_Fe
  i_Cu = i_BV(Ecorr-Erev(1), i0_Cu, alpha_Cu, z_Cu, temp);
  i_Fe = i_BV(Ecorr-Erev(2), i0_Fe, alpha_Fe, z_Fe, temp);
  func = i_Cu + i_Fe;
  end

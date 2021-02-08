%{
	Continuous stirred tank model for base metal extraction/recovery cell
    Using ode45
	To Add: Applied Voltage calc. Non-constant volume?
    Order of compounds: Inert Cu Sn Al Pb Fe Zn Ca Ni Ag Au Pd
%} 
clear variables

%universal constants
global R F;
R = 8.314; %J/(mol K)
F = 96485.3329; %C/mol
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
vol_lch = 200; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
m_PCB_total = 80; %kg Mass of crushed PCBs
r_particles = 0.001; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
tfinal = 259200; %s

%Weight fraction composition of PCB
%Inert Cu Sn Al Pb Fe 
wtfrac_PCB(1:6) = [0.783 0.151 1.80E-2 1.57E-2 1.16E-2 0.781E-2]; 
%Ag Au Pd
wtfrac_PCB(7:9) = [0.0130E-2 0.00580E-2 0.00286E-2];
%masses of each component:
m_PCB = m_PCB_total*wtfrac_PCB;
%Convert to volume and volume fraction
V_PCB = m_PCB./rho; %Element volumes in m^3
vfrac_PCB = V_PCB/sum(V_PCB);
V_PCB_total = sum(V_PCB);

packing_density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
fill_pct = 100*sum(V_PCB)/(0.001*packing_density*vol_lch);
if fill_pct > 75
    txt = ['Warning: PCB mass is ', num2str(fill_pct), '% of total lching vessel volume'];
    disp(txt);
end
%Surface area calculation for corrosion
SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
S_PCB_total = V_PCB_total*packing_density*SSA; %m2 Total surface area of spheres in lch.
S_PCB = vfrac_PCB*S_PCB_total; %array with exposed surface area of each component
S_PCB = S_PCB*100^2; %convert m^2 to cm^2

%Calculate number of crushed PCB particles
n_particles = sum(V_PCB)*3/(4*pi*r_particles^3);
%initial concentrations in mol/L
%Cell Concentrations (recovery)
Ci_Cu2_cell = 0;
Ci_Sn2_cell = 0.0;
Ci_Al3_cell = 0.0;
Ci_Pb2_cell = 0.0;
Ci_Fe2_cell = 0.5;
Ci_Fe3_cell = 0.1;
Ci_Ag_cell = 0.0;
Ci_Au_cell = 0.0;
Ci_Pd2_cell = 0.0;
Ci_H_cell = 0.5;
Ci_Cl_cell = 2*(Ci_Cu2_cell+Ci_Fe2_cell)+Ci_H_cell+3*Ci_Fe3_cell; %Calculation to ensure 
%electrolyte has net neutral charge
Ci_cell = [Ci_Cu2_cell Ci_Sn2_cell Ci_Al3_cell Ci_Pb2_cell Ci_Fe2_cell ...
    Ci_Fe3_cell Ci_Ag_cell Ci_Au_cell Ci_Pd2_cell Ci_H_cell Ci_Cl_cell];

%leching vessel concentrations (extraction)
Ci_Cu2_lch = 0;
Ci_Sn2_lch = 0;
Ci_Al3_lch = 0;
Ci_Pb2_lch = 0;
Ci_Fe2_lch = 0.5;
Ci_Fe3_lch = 0.1;
Ci_Ag_lch = 0;
Ci_Au_lch = 0;
Ci_Pd2_lch = 0;
Ci_H_lch = 0.5;
Ci_Cl_lch = 2*(Ci_Cu2_lch+Ci_Fe2_lch)+Ci_H_lch+3*Ci_Fe3_lch;
Ci_lch = [Ci_Cu2_lch Ci_Sn2_lch Ci_Al3_lch Ci_Pb2_lch Ci_Fe2_lch ... 
    Ci_Fe3_lch Ci_Ag_lch Ci_Au_lch Ci_Pd2_lch Ci_H_lch Ci_Cl_lch];

%initializing solution concentration and solid mass vector
Cm_i = [Ci_cell Ci_lch m_PCB];

%solve concentration profiles
tspan = [0 tfinal];
options = odeset('NonNegative',1:31);
balance_solver = @(t, Cm) ion_balance(t, Cm, temp, pres, vol_cell, vol_lch, Q, S_an, S_cat, V_app, n_particles, l, A_cell);
[t, Cm] = ode15s(balance_solver, tspan, Cm_i);
%{
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
    Erev_Cu_cell(j) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(Cm(j,1),eps)));
    Erev_Fe_cell(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(Cm(j,2),eps)/gamma_Fe3*max(Cm(j,3),eps));
    Erev_Sn_cell(j) = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma_Sn2*max(Cm(j,4),eps)));
    Erev_Ni_cell(j) = Eo_Ni - R*temp/(z_Ni*F)*log(1/(gamma_Ni2*max(Cm(j,5),eps)));
    Erev_Cu_bed(j) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*max(Cm(j,8),eps)));
    Erev_Fe_bed(j) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*max(Cm(j,9),eps)/gamma_Fe3*max(Cm(j,10),eps));
    Erev_Sn_bed(j) = Eo_Sn - R*temp/(z_Sn*F)*log(1/(gamma_Sn2*max(Cm(j,11),eps)));
    Erev_Ni_bed(j) = Eo_Ni - R*temp/(z_Ni*F)*log(1/(gamma_Ni2*max(Cm(j,12),eps)));
    
    %resistance calculation for IR drops
    kappa = 1000*(Cm(j,1)*lamda_Cu2 + Cm(j,2)*lamda_Fe2 + Cm(j,3)*lamda_Fe3 + Cm(j,4)*lamda_Sn2 + Cm(j,5)*lamda_Ni2 + Cm(j,6)*lamda_H + Cm(j,7)*lamda_Cl);
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
plot(t,Cm(:,1),t,Cm(:,8))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Copper')
subplot(3,2,2)
plot(t,Cm(:,2),t,Cm(:,9))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(II)')
subplot(3,2,3)
plot(t,Cm(:,3),t,Cm(:,10))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(III)')
subplot(3,2,4)
plot(t,Cm(:,4),t,Cm(:,11))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Tin')
subplot(3,2,5)
plot(t,Cm(:,5),t,Cm(:,12))
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
%}
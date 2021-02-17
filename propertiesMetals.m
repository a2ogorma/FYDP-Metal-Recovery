 %universal constants
global R F e;
R = 8.314; %J/(mol K)
F = 96485.3329; %C/mol
e = eps; %eps;realmin;
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
z_An = 4;
z = [z_Cu z_Sn z_Al z_Pb z_Fe1 z_Fe2 z_Ag z_Au z_Pd z_H z_An];
% exchange current densities
global i0
if solution == 1 %Cl-, base metal system
    active = [1 1 1 1 1 1 0 0 0 1 1];
else %S2O3, precious metal system
    active = [0 0 0 0 1 0 1 1 1 1 1];
end
i0_Cu = 5E-5; %A/cm^2
i0_Sn = 5E-6; %A/cm^2
i0_Al = 5E-6; %A/cm^2
i0_Pb = 5E-6; %A/cm^2
i0_Fe1 = 5E-6; %A/cm^2
i0_Fe2 = 5E-9; %A/cm^2
i0_Ag = 4e-9; %using a mix of a few from that reference
i0_Au = 2E-8; %A/cm2, value of dissolution including Na2S in dissolution reaction from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1134/S1023193506040021, converting A/cm2 to A/m2 
i0_Pd = 7.3e-9; %A/cm2, from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1007/s10800-007-9434-x - verify this link applies
i0_H = 1E-10; %need to identify suitable electrode to help limit hydrogen evolution. From volcano plot in class, likely copper, gold or silver electrodes
i0_An = 1E-12; %A/cm^2
i0 = active.*[i0_Cu i0_Sn i0_Al i0_Pb i0_Fe1 i0_Fe2 i0_Ag i0_Au i0_Pd i0_H i0_An];

global km
km = ones(1, 11)/(0.5*17367300); %mass transfer coefficient for limiting current (roughly results in limiting currents of 0.02)

% charge transfer coefficients
global alphas
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
alpha_An = 0.5;
alphas = [alpha_Cu alpha_Sn alpha_Al alpha_Pb alpha_Fe1 alpha_Fe2...
    alpha_Ag alpha_Au alpha_Pd alpha_H alpha_An];

% Standard half reaction potentials vs. SHE @ 298 K, 1 atm, 1 M https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
global Eo
Eo_Cu = 0.337; %V
Eo_Sn = -0.13; %V
Eo_Al = -1.66; %V
Eo_Pb = -0.126; %V
Eo_Fe1 = 0.77; %V
Eo_Fe2 = -0.44; %V
Eo_Ag = 0.060113; %V 0.7896; %V 
Eo_Au = 0.153; %V 1.83; %V 
Eo_Pd = 0.0862; %V 0.915; %V 
Eo_H = 0; %V
Eo_An = 1.23; %V
Eo = [Eo_Cu Eo_Sn Eo_Al Eo_Pb Eo_Fe1 Eo_Fe2 Eo_Ag Eo_Au Eo_Pd Eo_H Eo_An];

%mass transfer coefficients

%activity coefficients of ions in solution
global gamma aO2 aH2
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
gamma_S2O3 = 1;
gamma = [gamma_Cu2 gamma_Sn2 gamma_Al3 gamma_Pb2 gamma_Fe2 gamma_Fe3...
    gamma_Ag gamma_Au gamma_Pd2 gamma_H gamma_Cl gamma_S2O3];
aH2 = 0.0001; %atm
aO2 = 0.008/32.01; %Molar solubility in water at 0.21 atm partial pressue 

%Lamda infinity values NEED source for iron, tin, nickel rest are from ChE 331 notes
global lamda
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
lamda_S2O3 = 100; %S m^2/mol
lamda = [lamda_Cu2 lamda_Sn2 lamda_Al3 lamda_Pb2 lamda_Fe2 lamda_Fe3...
    lamda_Ag lamda_Au lamda_Pd lamda_H lamda_Cl lamda_S2O3];

%Density of inert material SUBJECT TO CHANGE (Based off range of epoxy)
global rho
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
global mw
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

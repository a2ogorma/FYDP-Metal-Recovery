%universal constants
global R F e;
R = 8.314; %J/(mol K)
F = 96485.3329; %C/mol
e = 1E-25; %eps;
%{
Ionic species order:
Base: Cu2+, Sn2+, Fe2+, Fe3+, Ag+, Au3+, Pd2+, H+, Cl-, (AuCl4)-,
Precious: Cu2+, Sn2+, Fe2+, Fe3+, [Ag(S2O3)2]3-, [Au(S2O3)2]3-,
    [Pd(S2O3)4]6-, H+, (S2O3)2-, (AuCl4)-
Solid mass order: Inert (If applicable), Cu, Sn, Fe, Ag, Au, Pd
%}
%{
Reactions
Cu2+ + 2e- <--> Cu(s) (1)
Sn2+ + 2e- <--> Sn(s) (2)
Fe3+ + e- <--> Fe2+ (3)
Fe2+ + 2e- <--> Fe(s) (4)
Ag+ + e- <--> Ag(s) (5B) OR [Ag(S2O3)2]3- + e- <--> Ag(s) + 2(S2O3)2- (5P)
AgCl(s) + e- <--> Ag(s) + Cl- (6)
Au3+ + 3e- <--> Au(s) (7B) OR [Au(S2O3)2]3- + e- <--> Au(s) + 2(S2O3)2- (7P)
[AuCl4]- + 3e- <--> Au(s) + 4Cl- (8)
Pd2+ + 2e- <--> Pd(s) (9B) OR [Pd(S2O3)4]6- + 2e- <--> Pd(s) + 4(S2O3)2- (9P)
H+ + e- <--> 0.5H2(g) (10)
H+ + 0.25O2(g) + e- <--> 0.5H2O(l) (11)
%}

%Electrochemical Constants
% # of electrons in rxn per mol of cation species consumed
global z;
if solution == 1 %Cl-, base metal system
    z = [2 2 1 2 1 1 3 3 2 1 1];
elseif solution == 2 %S2O3, precious metal system
    z = [2 2 1 2 1 1 1 3 2 1 1];
end

% exchange current densities, A/cm^2
if solution == 1 %Cl-, base metal system
    active = [1 0 1 1 0 0 0 0 0 1 1];
    i0_5 = 4e-9;
    i0_7 = 2E-8;
    i0_9 = 7.3E-9;
else %S2O3, precious metal system
    active = [1 1 1 1 1 1 1 0 1 1 1];
    i0_5 = 4e-9; %Aidan fill these out for the thiosulfate rxns pls
    i0_7 = 2E-8; %value of dissolution including Na2S in dissolution 
        %reaction from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1134/S1023193506040021, converting A/cm2 to A/m2 
    i0_9 = 7.3E-9;
end
i0_1 = 5E-5; %Cu
i0_2 = 5E-6; %Sn
i0_3 = 5E-6; %Fe3-> Fe2
i0_4 = 7E-6; %Fe2 -> Fe(s)
i0_6 = 4e-9; %AgCl
i0_8 = 2E-8; %AuCl 
i0_10 = 1E-10; % H2. need to identify suitable electrode to help limit hydrogen evolution. 
    %From volcano plot in class, likely copper, gold or silver electrodes
i0_11 = 1E-12; %OER
global i0
i0 = active.*[i0_1 i0_2 i0_3 i0_4 i0_5 i0_6 i0_7 i0_8 i0_9 i0_10 i0_11];

% charge transfer coefficients, m/s
global alphas
alphas = ones(1, 11)*0.5; %assume symmetric rxns 

% Standard half reaction potentials, V vs. SHE @ 298 K, 1 atm, 1 M conc.
    %https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
Eo_1 = 0.337;
Eo_2 = -0.13;
Eo_3 = 0.77;
Eo_4 = -0.44;
Eo_6 = 0.2223;
Eo_8 = 0.93;
Eo_10 = 0;
Eo_11  = 1.229;
if solution == 1 %Cl-, base metal system
    Eo_5 = 0.7996;
    Eo_7 = 1.52;
    Eo_9 = 0.915;
else %S2O3, precious metal system
    Eo_5 = 0.060113;
    Eo_7 = 0.153; 
    Eo_9 = 0.0862;
end
global Eo
Eo = [Eo_1 Eo_2 Eo_3 Eo_4 Eo_5 Eo_6 Eo_7 Eo_8 Eo_9 Eo_10 Eo_11];

%activity coefficients of ions in solution
global gamma aO2 aH2

gamma = ones(1,10);
aH2 = 0.0001; %atm
aO2 = 0.008/32.01; %Molar solubility in water at 0.21 atm partial pressue 

% ionic mass transfer coefficients
global km
km = ones(1, 10)*1.25E-3; %mass transfer coefficient for limiting current 
%based on km for Co3+ given in 331 notes

%Lamda infinity values NEED source for iron, tin, nickel rest are from ChE 331 notes
lamda_1 = 107.2; %Cu2
lamda_2 = 100; %Sn2
lamda_3 = 100; %Fe2
lamda_4 = 100; %Fe3
lamda_8 = 349.81; %H
lamda_10 = 100; %[AgCl4]
if solution == 1 %Cl-, base metal system
    lamda_5 = 100; %Ag
    lamda_6 = 100; %Au
    lamda_7 = 100; %Pd2
    lamda_9 = 76.35; %Cl-
else
    %thiosulfate complexes:
    lamda_5 = 100; %Ag
    lamda_6 = 100; %Au
    lamda_7 = 100; %Pd2
    lamda_9 = 100; %S2o3
end
global lamda
lamda = [lamda_1 lamda_2 lamda_3 lamda_4 lamda_5 lamda_6 lamda_7 lamda_8 ...
    lamda_9 lamda_10];

%Density of inert material SUBJECT TO CHANGE (Based off range of epoxy)
global rho
rho_inert = 1200; %kg/m^3
%metal densities in kg/m^3 
rho_Cu = 8960;
rho_Sn = 5750;
rho_Fe = 7847;
rho_Ag = 10490;
rho_Au = 19300;
rho_Pd = 11900;
rho = [rho_inert rho_Cu rho_Sn rho_Fe rho_Ag rho_Au rho_Pd];
%molecular weights
global mw
mw_Cu = 64.546;
mw_Sn = 118.71;
mw_Fe = 55.845;
mw_Ag = 107.8682;
mw_Au = 196.966;
mw_Pd = 106.42;
%index 1 empty to mark place of inert material in other arrays
mw = [0 mw_Cu mw_Sn mw_Fe mw_Ag mw_Au mw_Pd];

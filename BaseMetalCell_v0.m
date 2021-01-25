%{
	Continuous stirred tank model for base metal extraction/recovery cell
    Using ode45
	To Add: Applied Voltage calc. Non-constant volume?
%} 
clear variables

% system parameters
temp = 298; %K
pres = 1; % atm
vol_cell = 250; %L
Q = 5; % L/s (flowrate)
%Electrode areas
S_cat = 500; %cm^2
S_an = 500; %cm^2
%Applied Voltage (potentiostat)
V_app = 12; %V
%Extraction vessel parameters
vol_bed = 0.2; %m3 (Initial) volume of bed holding the particles assuming the bed is completly full.
r_particles = 0.001; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.

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
tspan = [0 10];
options = odeset('NonNegative', 1:10);
balance_solver = @(t, C) ion_balance(t, C, temp, pres, vol_cell, vol_bed, Q, S_an, S_cat, V_app, r_particles);
[t, C] = ode45(balance_solver, tspan, Ci, options);

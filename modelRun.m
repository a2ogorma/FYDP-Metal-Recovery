clear all
%% Base metal system parameters
solution = 1; %1 is Cl- base metal, 2 is S2O3 precious metal
propertiesMetals;

paramSetBase = struct;
paramSetBase.temp = 298; %K
paramSetBase.pres = 1; % atm
paramSetBase.Q = 0.001/60;%; % L/s (flowrate)
%cell dimension information
paramSetBase.length = 1; % m length of electrodes in flow direction x
paramSetBase.height = 0.5; % m height of electrodes
paramSetBase.spacing_x = 0.1; % m gap between end of electrode and vessel inlet/outlet
paramSetBase.spacing_y = 0.035; %m spacing between electrodes 
paramSetBase.n_epairs = 1; %number of anode cathode pairs
%paramSetBase.vol_cell = paramSetBase.n_epairs*2*paramSetBase.spacing_y*...
    paramSetBase.height*paramSetBase.length+2*paramSetBase.spacing_x;
paramSetBase.vol_cell = 0.04; % L
%Electrode areas
paramSetBase.S_cat = 2*40*6*28*pi*0.03302/2.54; %cm^2, area of 6 x 28 cm wire mesh
%paramSetBase.S_cat = paramSetBase.height*paramSetBase.length;
paramSetBase.S_an = 36; %cm^2
%paramSetBase.S_an = paramSetBase.S_cat;
%Cross sectional area of cell
paramSetBase.A_cell = 36; %cm^2
%paramSetBase.A_cell = paramSetBase.S_cat;
%L (Initial) volume of bed holding the particles assuming the bed is completly full.
paramSetBase.vol_lch = 0.09; 

paramSetBase.mode = 2; %1 - potentiostat, 2 - galvanostat
%Applied Voltage (potentiostat)
paramSetBase.V_app = 3; %V
%Applied Current to Cell (Galvanostat)
paramSetBase.I_app = 0.5; %A
paramSetBase.tfinal = 35*3600; %s


%Max current density for all rxns
paramSetBase.iL_default = 1; %A/cm^2
%fsolve options
paramSetBase.foptions = optimoptions(@fsolve, 'Display','off', ...
    'MaxFunctionEvaluations', 5000, 'Algorithm', 'trust-region-dogleg', 'StepTolerance', 1E-7);

%cathode calculations
thicc_cat = 0.01; %m, thickness of cathode
A_cat = paramSetBase.S_cat/2*10^-4;  %area of one side in m^2
s = A_cat^0.5; %assuming square shape
V_cat = A_cat*thicc_cat; %cathode volume in m^3
m_cat = rho(6)*V_cat; %Iron mass

%% Initial concentration specifications
initSetBase = struct;
%characteristics of solid PCB input
initSetBase.solidPCB.m_PCB_total = 0.07;%kg Mass of crushed PCBs
initSetBase.solidPCB.r_particles = 0.001;%m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
%Weight fraction composition of PCB
%Inert Cu Sn Fe 
initSetBase.solidPCB.wtfrac_PCB(1:4) = [0.783 0.151 1.80E-2 0.781E-2];  
%Ag Au Pd
initSetBase.solidPCB.wtfrac_PCB(5:7) = [0.0130E-2 0.00580E-2 0.00286E-2];
initSetBase.m_deposited = [0 0 m_cat 0 0 0]; 

%characteristics of starting solution
initSetBase.solution.type = solution;%1 is Cl- base metal, 2 is S2O3 precious metal
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSetBase.solution.Ci_Cu2_cell = 0.01;
initSetBase.solution.Ci_Sn2_cell = 0.0;
initSetBase.solution.Ci_Fe2_cell = 0.01;
initSetBase.solution.Ci_Fe3_cell = 0.5;
initSetBase.solution.Ci_Ag_cell = 0.00;
initSetBase.solution.Ci_Au3_cell = 0.0;
initSetBase.solution.Ci_Pd2_cell = 0.0;
initSetBase.solution.Ci_H_cell = 0.5;
%Calculation to ensure electrolyte has net neutral charge
initSetBase.solution.Ci_Cl_cell = 2*(initSetBase.solution.Ci_Cu2_cell+initSetBase.solution.Ci_Fe2_cell)+initSetBase.solution.Ci_H_cell+3*initSetBase.solution.Ci_Fe3_cell;
initSetBase.solution.Ci_AuCl4_cell = 0.0;
initSetBase.solution.Ci_cell = [initSetBase.solution.Ci_Cu2_cell initSetBase.solution.Ci_Sn2_cell initSetBase.solution.Ci_Fe2_cell ...
    initSetBase.solution.Ci_Fe3_cell initSetBase.solution.Ci_Ag_cell initSetBase.solution.Ci_Au3_cell ...
    initSetBase.solution.Ci_Pd2_cell initSetBase.solution.Ci_H_cell initSetBase.solution.Ci_Cl_cell initSetBase.solution.Ci_AuCl4_cell];

%leching vessel concentrations (extraction)
initSetBase.solution.Ci_Cu2_lch = initSetBase.solution.Ci_Cu2_cell;
initSetBase.solution.Ci_Sn2_lch = initSetBase.solution.Ci_Sn2_cell;
initSetBase.solution.Ci_Fe2_lch = initSetBase.solution.Ci_Fe2_cell;
initSetBase.solution.Ci_Fe3_lch = initSetBase.solution.Ci_Fe3_cell;
initSetBase.solution.Ci_Ag_lch = initSetBase.solution.Ci_Ag_cell;
initSetBase.solution.Ci_Au3_lch = initSetBase.solution.Ci_Au3_cell;
initSetBase.solution.Ci_Pd2_lch = initSetBase.solution.Ci_Pd2_cell;
initSetBase.solution.Ci_H_lch = initSetBase.solution.Ci_H_cell;
initSetBase.solution.Ci_Cl_lch = 2*(initSetBase.solution.Ci_Cu2_lch+initSetBase.solution.Ci_Fe2_lch)+initSetBase.solution.Ci_H_lch+3*initSetBase.solution.Ci_Fe3_lch;
initSetBase.solution.Ci_AuCl4_lch = initSetBase.solution.Ci_AuCl4_cell;
initSetBase.solution.Ci_lch = [initSetBase.solution.Ci_Cu2_lch initSetBase.solution.Ci_Sn2_lch initSetBase.solution.Ci_Fe2_lch ... 
    initSetBase.solution.Ci_Fe3_lch initSetBase.solution.Ci_Ag_lch initSetBase.solution.Ci_Au3_lch ...
    initSetBase.solution.Ci_Pd2_lch initSetBase.solution.Ci_H_lch initSetBase.solution.Ci_Cl_lch initSetBase.solution.Ci_AuCl4_lch];


%%
disp("Modelling Base Metal Extraction and Recovery");
resultsBase = metalER(initSetBase,paramSetBase);

%% Precious metal section (for connection)
solution = 2; %1 is Cl- base metal, 2 is S2O3 precious metal
propertiesMetals;
z = numel(resultsBase.t);
%characteristics of solid PCB input
initSetPrecious.solidPCB.m_PCB_total = resultsBase.PCB.massTotal(z);%kg Mass of crushed PCBs
initSetPrecious.solidPCB.r_particles = resultsBase.PCB.r_particles(z);%m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
initSetPrecious.m_deposited = [0 0 0 0.1 0 0];
%Weight fraction composition of PCB
%Inert Cu Sn Al Pb Fe 
initSetPrecious.solidPCB.wtfrac_PCB = resultsBase.PCB.wtfrac_PCB(z,:);

%characteristics of starting solution
initSetPrecious.solution.type = 2;%1 is Cl- base metal, 2 is S2O3 precious metal
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSetPrecious.solution.Ci_Cu2_cell = 0.01;
initSetPrecious.solution.Ci_Sn2_cell = 0.0;
initSetPrecious.solution.Ci_Fe2_cell = 0.1;
initSetPrecious.solution.Ci_Fe3_cell = 0.5;
initSetPrecious.solution.Ci_Ag_cell = 0.0;
initSetPrecious.solution.Ci_Au3_cell = 0.0;
initSetPrecious.solution.Ci_Pd2_cell = 0.0;
initSetPrecious.solution.Ci_H_cell = 0.5;
initSetPrecious.solution.Ci_S2O3_cell = 0.5;
initSetPrecious.solution.Ci_AuCl4_cell = 0;
initSetPrecious.solution.Ci_cell = [initSetPrecious.solution.Ci_Cu2_cell initSetPrecious.solution.Ci_Sn2_cell initSetPrecious.solution.Ci_Fe2_cell ...
    initSetPrecious.solution.Ci_Fe3_cell initSetPrecious.solution.Ci_Ag_cell initSetPrecious.solution.Ci_Au3_cell initSetPrecious.solution.Ci_Pd2_cell ...
    initSetPrecious.solution.Ci_H_cell initSetPrecious.solution.Ci_S2O3_cell initSetPrecious.solution.Ci_AuCl4_cell];

%leching vessel concentrations (extraction)
initSetPrecious.solution.Ci_Cu2_lch = 0;
initSetPrecious.solution.Ci_Sn2_lch = 0;
initSetPrecious.solution.Ci_Fe2_lch = 0.5;
initSetPrecious.solution.Ci_Fe3_lch = 0.1;
initSetPrecious.solution.Ci_Ag_lch = 0;
initSetPrecious.solution.Ci_Au3_lch = 0;
initSetPrecious.solution.Ci_Pd2_lch = 0;
initSetPrecious.solution.Ci_H_lch = 0.5;
initSetPrecious.solution.Ci_S2O3_lch = 0.5;
initSetPrecious.solution.Ci_AuCl4_lch = 0;
initSetPrecious.solution.Ci_lch = [initSetPrecious.solution.Ci_Cu2_lch initSetPrecious.solution.Ci_Sn2_lch initSetPrecious.solution.Ci_Fe2_lch ... 
    initSetPrecious.solution.Ci_Fe3_lch initSetPrecious.solution.Ci_Ag_lch initSetPrecious.solution.Ci_Au3_lch initSetPrecious.solution.Ci_Pd2_lch ...
    initSetPrecious.solution.Ci_H_lch initSetPrecious.solution.Ci_S2O3_lch initSetPrecious.solution.Ci_AuCl4_lch];

%% Precious metal system parameters
paramSetPrecious.temp = 298; %K
paramSetPrecious.pres = 1; % atm
paramSetPrecious.vol_cell = 0.04; %L
paramSetPrecious.Q = 0.001/60;%; % L/s (flowrate)
%Electrode areas
paramSetPrecious.S_cat = 250; %cm^2
paramSetPrecious.S_an = 36; %cm^2
%Cross sectional area of cell
paramSetPrecious.A_cell = 36; %cm^2
%Length b/w electrodes
paramSetPrecious.l = 3.5; %cm
%Extraction vessel parameters
paramSetPrecious.vol_lch = 0.09; %L (Initial) volume of bed holding the particles assuming the bed is completly full.

paramSetPrecious.mode = 1; %1 - potentiostat, 2 - galvanostat
%Applied Voltage (potentiostat)
paramSetPrecious.V_app = 3; %V
%Applied Current to Cell (Galvanostat)
paramSetPrecious.I_app = 0.05; %A
paramSetPrecious.tfinal = 10*3600; %s

%Max current density for all rxns
paramSetPrecious.iL_default = -1; %A/cm^2
%fsolve options
paramSetPrecious.foptions = optimoptions(@fsolve, 'Display','off', ...
    'MaxFunctionEvaluations', 5000, 'Algorithm', 'trust-region-dogleg', 'StepTolerance', 1E-7);

%%
%disp("Modelling Precious Metal Extraction and Recovery");
%resultsPrecious = metalER(initSetPrecious,paramSetPrecious);
%}
clear all

%% Precious metal section (for connection)
solution = 2; %1 is Cl- base metal, 2 is S2O3 precious metal
propertiesMetals;
%characteristics of solid PCB input
initSetPrecious.solidPCB.m_PCB_total = 5759;%kg Mass of crushed PCBs
initSetPrecious.solidPCB.r_particles = 0.0089;%m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
initSetPrecious.m_deposited = [eps 0 0 0 0 0];
%Weight fraction composition of PCB
%Inert Cu Sn Al Pb Fe 
wtfrac = [0.0845266104991119	0.2*0.815039475438435	0.2*0.0971570235621975	0.5*0.00210776765005768 0.000701689614615871 0.000313061520367081	0.000154371715215492];
initSetPrecious.solidPCB.wtfrac_PCB = wtfrac/sum(wtfrac);

%characteristics of starting solution
initSetPrecious.solution.type = 2;%1 is Cl- base metal, 2 is S2O3 precious metal
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSetPrecious.solution.Ci_Cu2_cell = 0.0;
initSetPrecious.solution.Ci_Sn2_cell = 0.0;
initSetPrecious.solution.Ci_Fe2_cell = 0.1;
initSetPrecious.solution.Ci_Fe3_cell = 0.5;
initSetPrecious.solution.Ci_Ag_cell = 0.0;
initSetPrecious.solution.Ci_Au3_cell = 0.00001;
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
solution = 2; %1 is Cl- base metal, 2 is S2O3 precious metal
propertiesMetals;
paramSetPrecious = struct;
paramSetPrecious.temp = 298; %K
paramSetPrecious.pres = 1; % atm
paramSetPrecious.Q = 47.8644;%0.0653/60;%; % L/s (flowrate)
%cell dimension information
paramSetPrecious.length = 4.7426; % m length of electrodes in flow direction x
paramSetPrecious.height = 4.1842; % m height of electrodes
paramSetPrecious.spacing_x = 0.1; % m gap between end of electrode and vessel inlet/outlet
paramSetPrecious.spacing_y = 0.1; %m spacing between electrodes 
paramSetPrecious.n_units = 45; %number of anode-cathode surface pairs
paramSetPrecious.vol_cell = (paramSetPrecious.n_units*paramSetPrecious.spacing_y*...
    paramSetPrecious.height*paramSetPrecious.length+2*paramSetPrecious.spacing_x)/1000; %L
%paramSetPrecious.vol_cell = 0.04; % L, volume of entire cell
%Electrode areas
%paramSetPrecious.S_cat = 36; %cm^2, area of 6 x 28 cm wire mesh
paramSetPrecious.S_cat = (paramSetPrecious.height/100)*(paramSetPrecious.length/100);
%paramSetPrecious.S_an = 36; %cm^2
paramSetPrecious.S_an = paramSetPrecious.S_cat;
%Cross sectional area of cell
%paramSetPrecious.A_cell = 36; %cm^2
paramSetPrecious.A_cell = paramSetPrecious.S_cat;
%L (Initial) volume of bed holding the particles assuming the bed is completly full.
paramSetPrecious.vol_lch = 6098.5; %L

paramSetPrecious.mode = 1; %1 - potentiostat, 2 - galvanostat
%Applied Voltage (potentiostat)
paramSetPrecious.V_app = 3; %V
%Applied Current to Cell (Galvanostat)
paramSetPrecious.I_app = 25;%36*0.01414; %A
paramSetPrecious.tfinal = 152180; %s

%Max current density for all rxns
paramSetPrecious.iL_default = 1; %A/cm^2
%fsolve options
paramSetPrecious.foptions = optimoptions(@fsolve, 'Display','off', ...
    'MaxFunctionEvaluations', 5000, 'Algorithm', 'trust-region-dogleg', 'StepTolerance', 1E-7);

%%
disp("Modelling Precious Metal Extraction and Recovery");
resultsPrecious = metalER(initSetPrecious,paramSetPrecious);
%}
metalViewStage2
clear all
%% Initial concentration specifications
initSet = struct;
%characteristics of solid PCB input
initSet.solidPCB.m_PCB_total = 80;%kg Mass of crushed PCBs
initSet.solidPCB.r_particles = 0.001;%m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
%Weight fraction composition of PCB
%Inert Cu Sn Al Pb Fe 
initSet.solidPCB.wtfrac_PCB(1:6) = [0.783 0.151 1.80E-2 1.57E-2 1.16E-2 0.781E-2];  
%Ag Au Pd
initSet.solidPCB.wtfrac_PCB(7:9) = [0.0130E-2 0.00580E-2 0.00286E-2];

%characteristics of starting solution
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSet.solution.Ci_Cu2_cell = 0;
initSet.solution.Ci_Sn2_cell = 0.0;
initSet.solution.Ci_Al3_cell = 0.0;
initSet.solution.Ci_Pb2_cell = 0.0;
initSet.solution.Ci_Fe2_cell = 0.5;
initSet.solution.Ci_Fe3_cell = 0.1;
initSet.solution.Ci_Ag_cell = 0.0;
initSet.solution.Ci_Au_cell = 0.0;
initSet.solution.Ci_Pd2_cell = 0.0;
initSet.solution.Ci_H_cell = 0.5;
initSet.solution.Ci_Cl_cell = 2*(initSet.solution.Ci_Cu2_cell+initSet.solution.Ci_Fe2_cell)+initSet.solution.Ci_H_cell+3*initSet.solution.Ci_Fe3_cell; %Calculation to ensure 
%electrolyte has net neutral charge
initSet.solution.Ci_cell = [initSet.solution.Ci_Cu2_cell initSet.solution.Ci_Sn2_cell initSet.solution.Ci_Al3_cell initSet.solution.Ci_Pb2_cell initSet.solution.Ci_Fe2_cell ...
    initSet.solution.Ci_Fe3_cell initSet.solution.Ci_Ag_cell initSet.solution.Ci_Au_cell initSet.solution.Ci_Pd2_cell initSet.solution.Ci_H_cell initSet.solution.Ci_Cl_cell];

%leching vessel concentrations (extraction)
initSet.solution.Ci_Cu2_lch = 0;
initSet.solution.Ci_Sn2_lch = 0;
initSet.solution.Ci_Al3_lch = 0;
initSet.solution.Ci_Pb2_lch = 0;
initSet.solution.Ci_Fe2_lch = 0.5;
initSet.solution.Ci_Fe3_lch = 0.1;
initSet.solution.Ci_Ag_lch = 0;
initSet.solution.Ci_Au_lch = 0;
initSet.solution.Ci_Pd2_lch = 0;
initSet.solution.Ci_H_lch = 0.5;
initSet.solution.Ci_Cl_lch = 2*(initSet.solution.Ci_Cu2_lch+initSet.solution.Ci_Fe2_lch)+initSet.solution.Ci_H_lch+3*initSet.solution.Ci_Fe3_lch;
initSet.solution.Ci_lch = [initSet.solution.Ci_Cu2_lch initSet.solution.Ci_Sn2_lch initSet.solution.Ci_Al3_lch initSet.solution.Ci_Pb2_lch initSet.solution.Ci_Fe2_lch ... 
    initSet.solution.Ci_Fe3_lch initSet.solution.Ci_Ag_lch initSet.solution.Ci_Au_lch initSet.solution.Ci_Pd2_lch initSet.solution.Ci_H_lch initSet.solution.Ci_Cl_lch];

%% system parameters
paramSet = struct;
paramSet.temp = 298; %K
paramSet.pres = 1; % atm
paramSet.vol_cell = 250; %L
paramSet.Q = 5;%; % L/s (flowrate)
%Electrode areas
paramSet.S_cat = 500; %cm^2
paramSet.S_an = 500; %cm^2
%Cross sectional area of cell
paramSet.A_cell = 500; %cm^2
%Length b/w electrodes
paramSet.l = 100; %cm
%Applied Voltage (potentiostat)
paramSet.V_app = 3; %V
%Extraction vessel parameters
paramSet.vol_lch = 200; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
paramSet.tfinal = 6*60*60; %s

results = BaseMetalCell(initSet,paramSet)
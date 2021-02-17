clear all
%% Initial concentration specifications
initSetBase = struct;
%characteristics of solid PCB input
initSetBase.solidPCB.m_PCB_total = 80;%kg Mass of crushed PCBs
initSetBase.solidPCB.r_particles = 0.001;%m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
%Weight fraction composition of PCB
%Inert Cu Sn Al Pb Fe 
initSetBase.solidPCB.wtfrac_PCB(1:6) = [0.783 0.151 1.80E-2 1.57E-2 1.16E-2 0.781E-2];  
%Ag Au Pd
initSetBase.solidPCB.wtfrac_PCB(7:9) = [0.0130E-2 0.00580E-2 0.00286E-2];
initSetBase.m_deposited = [0 0 0 0 0 0 0 0];

%characteristics of starting solution
initSetBase.solution.type = 1;%1 is Cl- base metal, 2 is S2O3 precious metal
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSetBase.solution.Ci_Cu2_cell = 0;
initSetBase.solution.Ci_Sn2_cell = 0.0;
initSetBase.solution.Ci_Al3_cell = 0.0;
initSetBase.solution.Ci_Pb2_cell = 0.0;
initSetBase.solution.Ci_Fe2_cell = 0.5;
initSetBase.solution.Ci_Fe3_cell = 0.1;
initSetBase.solution.Ci_Ag_cell = 0.0;
initSetBase.solution.Ci_Au_cell = 0.0;
initSetBase.solution.Ci_Pd2_cell = 0.0;
initSetBase.solution.Ci_H_cell = 0.5;
initSetBase.solution.Ci_Cl_cell = 2*(initSetBase.solution.Ci_Cu2_cell+initSetBase.solution.Ci_Fe2_cell)+initSetBase.solution.Ci_H_cell+3*initSetBase.solution.Ci_Fe3_cell; %Calculation to ensure 
initSetBase.solution.Ci_S2O3_cell = 0;
%electrolyte has net neutral charge
initSetBase.solution.Ci_cell = [initSetBase.solution.Ci_Cu2_cell initSetBase.solution.Ci_Sn2_cell initSetBase.solution.Ci_Al3_cell initSetBase.solution.Ci_Pb2_cell initSetBase.solution.Ci_Fe2_cell ...
    initSetBase.solution.Ci_Fe3_cell initSetBase.solution.Ci_Ag_cell initSetBase.solution.Ci_Au_cell initSetBase.solution.Ci_Pd2_cell initSetBase.solution.Ci_H_cell initSetBase.solution.Ci_Cl_cell initSetBase.solution.Ci_S2O3_cell];

%leching vessel concentrations (extraction)
initSetBase.solution.Ci_Cu2_lch = 0;
initSetBase.solution.Ci_Sn2_lch = 0;
initSetBase.solution.Ci_Al3_lch = 0;
initSetBase.solution.Ci_Pb2_lch = 0;
initSetBase.solution.Ci_Fe2_lch = 0.5;
initSetBase.solution.Ci_Fe3_lch = 0.1;
initSetBase.solution.Ci_Ag_lch = 0;
initSetBase.solution.Ci_Au_lch = 0;
initSetBase.solution.Ci_Pd2_lch = 0;
initSetBase.solution.Ci_H_lch = 0.5;
initSetBase.solution.Ci_Cl_lch = 2*(initSetBase.solution.Ci_Cu2_lch+initSetBase.solution.Ci_Fe2_lch)+initSetBase.solution.Ci_H_lch+3*initSetBase.solution.Ci_Fe3_lch;
initSetBase.solution.Ci_S2O3_lch = 0;
initSetBase.solution.Ci_lch = [initSetBase.solution.Ci_Cu2_lch initSetBase.solution.Ci_Sn2_lch initSetBase.solution.Ci_Al3_lch initSetBase.solution.Ci_Pb2_lch initSetBase.solution.Ci_Fe2_lch ... 
    initSetBase.solution.Ci_Fe3_lch initSetBase.solution.Ci_Ag_lch initSetBase.solution.Ci_Au_lch initSetBase.solution.Ci_Pd2_lch initSetBase.solution.Ci_H_lch initSetBase.solution.Ci_Cl_lch initSetBase.solution.Ci_S2O3_lch];

%% Base metal system parameters
paramSetBase = struct;
paramSetBase.temp = 298; %K
paramSetBase.pres = 1; % atm
paramSetBase.vol_cell = 250; %L
paramSetBase.Q = 5;%; % L/s (flowrate)
%Electrode areas
paramSetBase.S_cat = 500; %cm^2
paramSetBase.S_an = 500; %cm^2
%Cross sectional area of cell
paramSetBase.A_cell = 500; %cm^2
%Length b/w electrodes
paramSetBase.l = 100; %cm
%Applied Voltage (potentiostat)
paramSetBase.V_app = 3; %V
%Extraction vessel parameters
paramSetBase.vol_lch = 200; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
paramSetBase.tfinal = 10000; %s
%%
resultsBase = BaseMetalCell(initSetBase,paramSetBase)

%% Precious metal section (for connection)
%{ 
%characteristics of solid PCB input
initSetPrecious.solidPCB.m_PCB_total = resultsBase.PCB.massTotal(end);%kg Mass of crushed PCBs
initSetPrecious.solidPCB.r_particles = resultsBase.PCB.r_particles(end);%m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
%Weight fraction composition of PCB
%Inert Cu Sn Al Pb Fe 
initSetPrecious.solidPCB.wtfrac_PCB = resultsBase.PCB.wtfrac_PCB(end,:);

%characteristics of starting solution
initSetPrecious.solution.type = 2;%1 is Cl- base metal, 2 is S2O3 precious metal
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSetPrecious.solution.Ci_Cu2_cell = 0;
initSetPrecious.solution.Ci_Sn2_cell = 0.0;
initSetPrecious.solution.Ci_Al3_cell = 0.0;
initSetPrecious.solution.Ci_Pb2_cell = 0.0;
initSetPrecious.solution.Ci_Fe2_cell = 0.5;
initSetPrecious.solution.Ci_Fe3_cell = 0.1;
initSetPrecious.solution.Ci_Ag_cell = 0.0;
initSetPrecious.solution.Ci_Au_cell = 0.0;
initSetPrecious.solution.Ci_Pd2_cell = 0.0;
initSetPrecious.solution.Ci_H_cell = 0.5;
initSetPrecious.solution.Ci_Cl_cell = 0;
initSetPrecious.solution.Ci_S2O3_cell = 0.5;
initSetPrecious.solution.Ci_cell = [initSetPrecious.solution.Ci_Cu2_cell initSetPrecious.solution.Ci_Sn2_cell initSetPrecious.solution.Ci_Al3_cell initSetPrecious.solution.Ci_Pb2_cell initSetPrecious.solution.Ci_Fe2_cell ...
    initSetPrecious.solution.Ci_Fe3_cell initSetPrecious.solution.Ci_Ag_cell initSetPrecious.solution.Ci_Au_cell initSetPrecious.solution.Ci_Pd2_cell initSetPrecious.solution.Ci_H_cell initSetPrecious.solution.Ci_Cl_cell initSetPrecious.solution.Ci_S2O3_cell];

%leching vessel concentrations (extraction)
initSetPrecious.solution.Ci_Cu2_lch = 0;
initSetPrecious.solution.Ci_Sn2_lch = 0;
initSetPrecious.solution.Ci_Al3_lch = 0;
initSetPrecious.solution.Ci_Pb2_lch = 0;
initSetPrecious.solution.Ci_Fe2_lch = 0.5;
initSetPrecious.solution.Ci_Fe3_lch = 0.1;
initSetPrecious.solution.Ci_Ag_lch = 0;
initSetPrecious.solution.Ci_Au_lch = 0;
initSetPrecious.solution.Ci_Pd2_lch = 0;
initSetPrecious.solution.Ci_H_lch = 0.5;
initSetPrecious.solution.Ci_Cl_lch = 0;
initSetPrecious.solution.Ci_S2O3_lch = 0.5;
initSetPrecious.solution.Ci_lch = [initSetPrecious.solution.Ci_Cu2_lch initSetPrecious.solution.Ci_Sn2_lch initSetPrecious.solution.Ci_Al3_lch initSetPrecious.solution.Ci_Pb2_lch initSetPrecious.solution.Ci_Fe2_lch ... 
    initSetPrecious.solution.Ci_Fe3_lch initSetPrecious.solution.Ci_Ag_lch initSetPrecious.solution.Ci_Au_lch initSetPrecious.solution.Ci_Pd2_lch initSetPrecious.solution.Ci_H_lch initSetPrecious.solution.Ci_Cl_lch initSetPrecious.solution.Ci_S2O3_lch];

%% Precious metal system parameters
paramSetPrecious.temp = 298; %K
paramSetPrecious.pres = 1; % atm
paramSetPrecious.vol_cell = 250; %L
paramSetPrecious.Q = 5;%; % L/s (flowrate)
%Electrode areas
paramSetPrecious.S_cat = 500; %cm^2
paramSetPrecious.S_an = 500; %cm^2
%Cross sectional area of cell
paramSetPrecious.A_cell = 500; %cm^2
%Length b/w electrodes
paramSetPrecious.l = 100; %cm
%Applied Voltage (potentiostat)
paramSetPrecious.V_app = 7; %V
%Extraction vessel parameters
paramSetPrecious.vol_lch = 200; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
paramSetPrecious.tfinal = 6*60*60; %s
%%
resultsPrecious = BaseMetalCell(initSetPrecious,paramSetPrecious)
%}
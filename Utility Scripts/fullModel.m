clear all
%% Preprocessing %%
resultsPreprocessing.wtFracIn = [0.753622 0.1936 0.0231 0.0294 167E-6 74.33E-6 36.67E-6]; %Inert Cu Sn Fe Ag Au Pd
resultsPreprocessing.Throughput = 400000; %kg/yr, raw PCB feed
CF = 0.91; %capacity factor, operation hours per year
resultsPreprocessing.CF = CF;
resultsPreprocessing.workingFactor = 0.25; %percentage of annual hours the preprocessing system works
resultsPreprocessing.massFlowrate = resultsPreprocessing.Throughput/(resultsPreprocessing.CF*resultsPreprocessing.workingFactor*8760*3600); %kg/s
%grinder
resultsPreprocessing.d_particles = 0.002; %m, size coming out of grinder
resultsPreprocessing.grinder.power = 0.008*resultsPreprocessing.massFlowrate/resultsPreprocessing.d_particles;%kW
loss = 0.001; %fraction of material lost
resultsPreprocessing.grinder.output = (1-loss)*resultsPreprocessing.Throughput; %kg/yr
g_out = resultsPreprocessing.grinder.output*resultsPreprocessing.wtFracIn; %kg/yr, partial throughputs
%ESP
resultsPreprocessing.ESP.power = 31; %kW
nonmetalSeparationEff = 0.98; %fraction of inert material separated from the system
metalRecoveryEff = 1; % percentage of metals recovered from ESP
ESP_out = g_out.*[(1-nonmetalSeparationEff) 1 1 1 1 1 1]; %kg/yr, partial throughputs
resultsPreprocessing.ESP.wtFracOut = ESP_out/sum(ESP_out); %weight fraction.
resultsPreprocessing.ESP.mainOutput = sum(ESP_out); %kg/yr
resultsPreprocessing.ESP.wasteOutput = resultsPreprocessing.grinder.output - resultsPreprocessing.ESP.mainOutput; %kg/yr (non-metals)
%Drum
resultsPreprocessing.drum.power = 0.02; %kW
ferrousSeparationEff = 0.95;
drum_out = ESP_out; %kg/yr, partial throughputs
drum_out(4) = ESP_out(4)*(1-ferrousSeparationEff); %kg/yr, partial throughputs
resultsPreprocessing.drum.wtFracOut = drum_out/sum(drum_out); %weight frac out of drum
resultsPreprocessing.drum.mainOutput = sum(drum_out);%kg/yr
resultsPreprocessing.drum.ferroOutput = resultsPreprocessing.ESP.mainOutput - resultsPreprocessing.drum.mainOutput; %kg/yr
%final calcs
resultsPreprocessing.productionRate = resultsPreprocessing.drum.mainOutput; %kg/yr
resultsPreprocessing.wtFracOut = resultsPreprocessing.drum.wtFracOut;
ModelResults.resultsPreprocessing = resultsPreprocessing;
%% Base Metal Recovery %%

%Base metal system parameters
solution = 1; %1 is Cl- base metal, 2 is S2O3 precious metal
propertiesMetals;
initSetBase = struct;
%characteristics of solid PCB output from preprocessing
initSetBase.solidPCB.wtfrac_PCB = resultsPreprocessing.wtFracOut;
%Cycle time for Base and Metal recovery operations
tfinal = 30*3600;
%PCB mass loaded per cycle
initSetBase.solidPCB.m_PCB_total = resultsPreprocessing.productionRate/(8760*CF)*tfinal/3600; %kg/yr/(hr/yr)*hr/cycle = kg/cycle
global rho
V_PCB_total = sum(initSetBase.solidPCB.m_PCB_total.*initSetBase.solidPCB.wtfrac_PCB./rho)*1000;%L
%Particle radius, m
initSetBase.solidPCB.r_particles = resultsPreprocessing.d_particles/2; 

%characteristics of starting solution
initSetBase.solution.type = solution;%1 is Cl- base metal, 2 is S2O3 precious metal
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSetBase.solution.Ci_Cu2_cell = 0.0001;
initSetBase.solution.Ci_Sn2_cell = 0.0001;
initSetBase.solution.Ci_Fe2_cell = 0.001;
initSetBase.solution.Ci_Fe3_cell = 1;
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
%Base metal setup
paramSetBase = struct;
paramSetBase.temp = 298; %K
paramSetBase.pres = 1; % atm
paramSetBase.Q = 4;% L/s (flowrate)
%cell dimension information
paramSetBase.length = 6; % m length of electrodes in flow direction x
paramSetBase.height = 1; % m height of electrodes    paramSetBase.spacing_x = 0.1; % m gap between end of electrode and vessel inlet/outlet
paramSetBase.spacing_y = 0.045; %m spacing between electrodes 
paramSetBase.spacing_x = 0.05; %m spacing between end of cell and electrodes
paramSetBase.n_units = 5; %number of anode-cathode surface pairs
paramSetBase.vol_cell = (paramSetBase.n_units*paramSetBase.spacing_y*...
    paramSetBase.height*(paramSetBase.length+2*paramSetBase.spacing_x))*1000; %L Volume of electrolyte in cell
%Electrode areas (one side), cm^2
paramSetBase.S_cat = (paramSetBase.height*100)*(paramSetBase.length*100);
paramSetBase.S_an = paramSetBase.S_cat;
%Cross sectional area of cell
paramSetBase.A_cell = paramSetBase.S_cat;
%L (Initial) volume of bed holding the particles assuming the bed is half
%full
vol_bed = (V_PCB_total/0.6/0.5);
paramSetBase.vol_lch = vol_bed-V_PCB_total; %L, volume of electrolyte in bed

paramSetBase.mode = 1; %1 - potentiostat, 2 - galvanostat
%Applied Voltage (potentiostat)
paramSetBase.V_app = 3; %V
%Applied Current to Cell (Galvanostat)
paramSetBase.I_app = 36*0.01414; %A
%Processing time
paramSetBase.tfinal = tfinal; %s

%Max current density for all rxns
paramSetBase.iL_default = 1; %A*m/dm^3
%fsolve options
paramSetBase.foptions = optimoptions(@fsolve, 'Display','off', ...
    'MaxFunctionEvaluations', 5000, 'Algorithm', 'trust-region-dogleg', 'StepTolerance', 1E-7);

%cathode calculations
thicc_cat = 0.01; %m, thickness of cathode
A_cat = paramSetBase.S_cat/2*10^-4;  %area of one side in m^2
s = A_cat^0.5; %assuming square shape
V_cat = A_cat*thicc_cat; %cathode volume in m^3
m_cat = rho(6)*V_cat; %Iron mass
initSetBase.m_deposited = [0 0 m_cat 0 0 0]; %Cu Sn Fe Ag Au Pd


disp("Modelling Base Metal Extraction and Recovery");
resultsBase = metalER(initSetBase,paramSetBase);
%Post calculations for impact metrics
%Practical additions here that dont affect the model
resultsBase.practical.pump.flow = resultsBase.init.paramSet.Q; %flow rate in system
resultsBase.practical.pump.head = 3; %reasonable assumption value, m 
resultsBase.practical.pump.specGravity = 1;
resultsBase.practical.pump.shaftPower = 9.81*resultsBase.practical.pump.specGravity*resultsBase.practical.pump.flow*1000*resultsBase.practical.pump.head/1000;
resultsBase.practical.pump.eff = 0.5;
resultsBase.practical.pump.BHP = resultsBase.practical.pump.shaftPower/resultsBase.practical.pump.eff;
%final stuff
resultsBase.productionRateIn = resultsBase.init.initSet.solidPCB.m_PCB_total/resultsBase.init.paramSet.tfinal*(3600*8760*CF); %kg/s*s/yr = kg/yr
resultsBase.numberUnits = ceil(resultsPreprocessing.productionRate/resultsBase.productionRateIn);
resultsBase.productionRateOut = resultsBase.numberUnits*resultsBase.PCB.massTotal(end)/resultsBase.init.paramSet.tfinal*(3600*8760*CF); %kg/s*s/yr = kg/yr
ModelResults.resultsBase = resultsBase;

%% Precious Metal Recovery %%
solution = 2; %1 is Cl- base metal, 2 is S2O3 precious metal
propertiesMetals;
initSetPrecious = struct;
%characteristics of solid PCB input
initSetPrecious.solidPCB.r_particles = resultsBase.PCB.r_particles(end);
m_PCB_p = resultsBase.PCB.massRem(end,:);
initSetPrecious.solidPCB.m_PCB_total = sum(m_PCB_p);
initSetPrecious.solidPCB.wtfrac_PCB = m_PCB_p/sum(m_PCB_p);
initSetPrecious.m_deposited = [eps 0 0 0 0 0]; %Given a small mass of iron to prevent error
V_PCB_total = sum(initSetPrecious.solidPCB.m_PCB_total.*initSetPrecious.solidPCB.wtfrac_PCB./rho)*1000;%L

%characteristics of starting solution
initSetPrecious.solution.type = solution;%1 is Cl- base metal, 2 is S2O3 precious metal
%initial concentrations in mol/L
%Cell Concentrations (recovery)
initSetPrecious.solution.Ci_Cu2_cell = 0.0;
initSetPrecious.solution.Ci_Sn2_cell = 0.0;
initSetPrecious.solution.Ci_Fe2_cell = 0.1;
initSetPrecious.solution.Ci_Fe3_cell = 0.5;
initSetPrecious.solution.Ci_Ag_cell = 0.0;
initSetPrecious.solution.Ci_Au3_cell = 0.0;
initSetPrecious.solution.Ci_Pd2_cell = 0.0;
initSetPrecious.solution.Ci_H_cell = 1E-10;
initSetPrecious.solution.Ci_S2O3_cell = 0.5;
initSetPrecious.solution.Ci_AuCl4_cell = 0;
initSetPrecious.solution.Ci_cell = [initSetPrecious.solution.Ci_Cu2_cell initSetPrecious.solution.Ci_Sn2_cell initSetPrecious.solution.Ci_Fe2_cell ...
    initSetPrecious.solution.Ci_Fe3_cell initSetPrecious.solution.Ci_Ag_cell initSetPrecious.solution.Ci_Au3_cell initSetPrecious.solution.Ci_Pd2_cell ...
    initSetPrecious.solution.Ci_H_cell initSetPrecious.solution.Ci_S2O3_cell initSetPrecious.solution.Ci_AuCl4_cell];

%leching vessel concentrations (extraction)
initSetPrecious.solution.Ci_Cu2_lch = initSetPrecious.solution.Ci_Cu2_cell;
initSetPrecious.solution.Ci_Sn2_lch = initSetPrecious.solution.Ci_Sn2_cell;
initSetPrecious.solution.Ci_Fe2_lch = initSetPrecious.solution.Ci_Fe2_cell;
initSetPrecious.solution.Ci_Fe3_lch = initSetPrecious.solution.Ci_Fe3_cell;
initSetPrecious.solution.Ci_Ag_lch = initSetPrecious.solution.Ci_Ag_cell;
initSetPrecious.solution.Ci_Au3_lch = initSetPrecious.solution.Ci_Au3_cell;
initSetPrecious.solution.Ci_Pd2_lch = initSetPrecious.solution.Ci_Pd2_cell;
initSetPrecious.solution.Ci_H_lch = initSetPrecious.solution.Ci_H_cell;
initSetPrecious.solution.Ci_S2O3_lch = initSetPrecious.solution.Ci_S2O3_cell;
initSetPrecious.solution.Ci_AuCl4_lch = initSetPrecious.solution.Ci_AuCl4_cell;
initSetPrecious.solution.Ci_lch = [initSetPrecious.solution.Ci_Cu2_lch initSetPrecious.solution.Ci_Sn2_lch initSetPrecious.solution.Ci_Fe2_lch ... 
    initSetPrecious.solution.Ci_Fe3_lch initSetPrecious.solution.Ci_Ag_lch initSetPrecious.solution.Ci_Au3_lch initSetPrecious.solution.Ci_Pd2_lch ...
    initSetPrecious.solution.Ci_H_lch initSetPrecious.solution.Ci_S2O3_lch initSetPrecious.solution.Ci_AuCl4_lch];
paramSetPrecious = struct;
paramSetPrecious.temp = 298; %K
paramSetPrecious.pres = 1; % atm
paramSetPrecious.Q = 2;% L/s (flowrate)
%cell dimension information
paramSetPrecious.length = 5; % m length of electrodes in flow direction x
paramSetPrecious.height = 2; % m height of electrodes
paramSetPrecious.spacing_x = 0.1; % m gap between end of electrode and vessel inlet/outlet
paramSetPrecious.spacing_y = 0.045; %m spacing between electrodes 
paramSetPrecious.n_units = 5; %number of anode-cathode surface pairs
paramSetPrecious.vol_cell = (paramSetPrecious.n_units*paramSetPrecious.spacing_y*...
    paramSetPrecious.height*(paramSetPrecious.length+2*paramSetPrecious.spacing_x))*1000; %L
%Electrode areas, cm^2
paramSetPrecious.S_cat = (paramSetPrecious.height*100)*(paramSetPrecious.length*100);
paramSetPrecious.S_an = paramSetPrecious.S_cat;
%Cross sectional area of cell
paramSetPrecious.A_cell = paramSetPrecious.S_cat;
%L (Initial) volume of bed holding the particles assuming the bed is half
%full
vol_bed = (V_PCB_total/0.6/0.5);
paramSetPrecious.vol_lch = vol_bed-V_PCB_total; %L, volume of electrolyte in bed


paramSetPrecious.mode = 1; %1 - potentiostat, 2 - galvanostat
%Applied Voltage (potentiostat)
paramSetPrecious.V_app = 3; %V
%Applied Current to Cell (Galvanostat)
paramSetPrecious.I_app = 25;%36*0.01414; %A
paramSetPrecious.tfinal = tfinal; %s

%Max current density for all rxns
paramSetPrecious.iL_default = 1; %A/cm^2
%fsolve options
paramSetPrecious.foptions = optimoptions(@fsolve, 'Display','off', ...
    'MaxFunctionEvaluations', 5000, 'Algorithm', 'trust-region-dogleg', 'StepTolerance', 1E-7);

%run ER model:
disp("Modelling Precious Metal Extraction and Recovery");
resultsPrecious = metalER(initSetPrecious,paramSetPrecious);
%Practical additions here that dont affect the model
resultsPrecious.practical.pump.flow = resultsPrecious.init.paramSet.Q; %flow rate in system
resultsPrecious.practical.pump.head = 3; %reasonable assumption value
resultsPrecious.practical.pump.specGravity = 1;
resultsPrecious.practical.pump.shaftPower = 9.81*resultsPrecious.practical.pump.specGravity*resultsPrecious.practical.pump.flow*1000*resultsPrecious.practical.pump.head/1000;
resultsPrecious.practical.pump.eff = 0.5;
resultsPrecious.practical.pump.BHP = resultsPrecious.practical.pump.shaftPower/resultsPrecious.practical.pump.eff;
%final stuff
resultsPrecious.productionRateIn = resultsPrecious.init.initSet.solidPCB.m_PCB_total/resultsPrecious.init.paramSet.tfinal*(3600*8760*CF); %kg/s*s/yr = kg/yr
resultsPrecious.numberUnits = ceil(resultsBase.productionRateOut/resultsPrecious.productionRateIn);
resultsPrecious.productionRateOut = resultsPrecious.numberUnits*resultsPrecious.PCB.massTotal(end)/resultsPrecious.init.paramSet.tfinal*(3600*8760*CF); %kg/s*s/yr = kg/yr
ModelResults.resultsPrecious = resultsPrecious;

%% Impacts %%
[resultsEnvironmental, resultsEconomic] = impactMetrics(resultsPreprocessing, resultsBase, resultsPrecious);
ModelResults.resultsEnvironmental = resultsEnvironmental;
ModelResults.resultsEconomic = resultsEconomic;
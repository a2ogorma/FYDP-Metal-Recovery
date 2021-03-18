resultsPreprocessing.wtFracIn = [0.783 0.151 1.80E-2 0.781E-2 0.0130E-2 0.00580E-2 0.00286E-2]; %Inert Cu Sn Fe Ag Au Pd
resultsPreprocessing.Throughput = 100000; %kg
resultsPreprocessing.CF = 0.91; %capacity factor, operation hours per year
resultsPreprocessing.workingFactor = 0.25; %percentage of annual hours the preprocessing system works
resultsPreprocessing.massFlowrate = resultsPreprocessing.Throughput/(resultsPreprocessing.CF*resultsPreprocessing.workingFactor*8760*3600); %kg/s
%grinder
resultsPreprocessing.d_particles = 0.001; %m, size coming out of grinder
resultsPreprocessing.grinder.power = 0.008*resultsPreprocessing.massFlowrate/resultsPreprocessing.d_particles;
loss = 0.001; %amount of material lost
resultsPreprocessing.grinder.output = (1-loss)*resultsPreprocessing.massFlowrate;
%ESP
resultsPreprocessing.ESP.power = 31; %kW
nonmetalSeperationEff = 0.9972; %percentage of inert material separated from the system
metalRecoveryEff = 1; % percentage of metals recovered from ESP
resultsPreprocessing.ESP.wtFracOut = [resultsPreprocessing.wtFracIn(1)*(1-nonmetalSeperationEff) metalRecoveryEff.*resultsPreprocessing.wtFracIn(2:7)]; %portion of original stream
resultsPreprocessing.ESP.mainOutput = sum(resultsPreprocessing.Throughput*resultsPreprocessing.ESP.wtFracOut);
resultsPreprocessing.ESP.wasteOutput = resultsPreprocessing.Throughput - resultsPreprocessing.ESP.mainOutput;
%Drum
resultsPreprocessing.drum.power = 0.02; %kW
ferrousSeperationEff = 0.95;
resultsPreprocessing.drum.wtFracOut = resultsPreprocessing.ESP.wtFracOut; %setting up for next line
resultsPreprocessing.drum.wtFracOut(4) = resultsPreprocessing.drum.wtFracOut(4)*(1-ferrousSeperationEff);
resultsPreprocessing.drum.mainOutput = sum(resultsPreprocessing.Throughput*resultsPreprocessing.drum.wtFracOut);
resultsPreprocessing.drum.wasteOutput = resultsPreprocessing.ESP.mainOutput - resultsPreprocessing.drum.mainOutput;
%final calcs
resultsPreprocessing.productionRate = resultsPreprocessing.drum.wasteOutput/(resultsPreprocessing.CF*resultsPreprocessing.workingFactor*8760*3600); %in kg per second
resultsPreprocessing.wtFracOut = resultsPreprocessing.drum.wtFracOut/sum(resultsPreprocessing.drum.wtFracOut)
%% preprocessing sizing %%

ModelResults.resultsPreprocessing = resultsPreprocessing;
%% Base metals %%
initSetBase.solidPCB.r_particles = resultsPreprocessing.d_particles/2;
%Practical additions here that done affect the model
resultsBase.practical.pump.flow = resultsBase.init.paramSet.Q; %flow rate in system
resultsBase.practical.pump.head = 3; %reasonable assumption value
resultsBase.practical.pump.specGravity = 1;
resultsBase.practical.pump.shaftPower = 9.81*resultsBase.practical.pump.specGravity*resultsBase.practical.pump.flow*1000*resultsBase.practical.pump.head/1000;
resultsBase.practical.pump.eff = 0.5;
resultsBase.practical.pump.BHP = resultsBase.practical.pump.shaftPower/resultsBase.practical.pump.eff;
%final stuff
resultsBase.productionRateIn = resultsBase.init.initSet.solidPCB.m_PCB_total/resultsBase.init.paramSet.tfinal;
resultsBase.productionRateOut = resultsBase.PCB.massTotal(end)/resultsBase.init.paramSet.tfinal;
resultsBase.numberUnits = resultsPreprocessing.productionRate/resultsBase.productionRateIn;
ModelResults.resultsBase = resultsBase;
%% Precious metals %%
%Practical additions here that done affect the model
resultsPrecious.practical.pump.flow = resultsPrecious.init.paramSet.Q; %flow rate in system
resultsPrecious.practical.pump.head = 3; %reasonable assumption value
resultsPrecious.practical.pump.specGravity = 1;
resultsPrecious.practical.pump.shaftPower = 9.81*resultsPrecious.practical.pump.specGravity*resultsPrecious.practical.pump.flow*1000*resultsPrecious.practical.pump.head/1000;
resultsPrecious.practical.pump.eff = 0.5;
resultsPrecious.practical.pump.BHP = resultsPrecious.practical.pump.shaftPower/resultsPrecious.practical.pump.eff;
%final stuff
resultsPrecious.productionRateIn = resultsPrecious.init.initSet.solidPCB.m_PCB_total/resultsPrecious.init.paramSet.tfinal;
resultsPrecious.numberUnits = resultsBase.productionRateOut/resultsPrecious.productionRateIn;
resultsPrecious.productionRateOut = resultsPrecious.PCB.massTotal(end)/resultsPrecious.init.paramSet.tfinal;
ModelResults.resultsPrecious = resultsPrecious;
%% Impacts %%
[resultsEnvironmental, resultsEconomic] = impactMetrics(resultsPreprocessing, resultsBase, resultsPrecious);
ModelResults.resultsEnvironmental = resultsEnvironmental;
ModelResults.resultsEconomic = resultsEconomic;
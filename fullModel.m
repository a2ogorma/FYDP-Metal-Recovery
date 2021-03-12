Throughput = 100000; %kg
CF = 0.91; %capacity factor, operation hours per year
%% preprocessing sizing %%
resultsPreprocessing.productionRate = productionRate; %in kg per second
ModelResults.resultsPreprocessing = resultsPreprocessing;
%% Base metals %%
resultsBase.productionRate = resultsBase.init.initSet.solidPCB.m_PCB_total/resultsBase.init.paramSet.tfinal;
resultsBase.numberUnits = resultsPreprocessing.productionRate/resultsBase.productionRate;
ModelResults.resultsBase = resultsBase;
%% Precious metals %%
resultsPrecious.productionRate = resultsPrecious.init.initSet.solidPCB.m_PCB_total/resultsPrecious.init.paramSet.tfinal;
resultsPrecious.numberUnits = resultsBase.productionRate/resultsPrecious.productionRate;
ModelResults.resultsPrecious = resultsPrecious;
%% Impacts %%
resultsEnvironmental = environmentalImpact(resultsPreprocessing, resultsBase, resultsPrecious);
ModelResults.resultsEnvironmental = resultsEnvironmental;

resultsEconomics = Economics(resultsPreprocessing, resultsPrecious, resultsBase);
ModelResults.resultsEconomics = resultsEconomics;
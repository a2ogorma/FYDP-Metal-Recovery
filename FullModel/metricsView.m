clear all
%% Open Files
sim_files = dir(fullfile('FullModel', '*.mat'));
ModelResults = open(fullfile('FullModel',sim_files(1).name)).ModelResults;
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('FullModel',sim_files(k).name)).ModelResults;
end

%% Extract Data
for n = 1:1:length(ModelResults)
    paramSetBase(n) = ModelResults(n).resultsBase.init.paramSet;
    initSetBase(n) = ModelResults(n).resultsBase.init.initSet;
    paramSetPrecious(n) = ModelResults(n).resultsPrecious.init.paramSet;
    initSetPrecious(n) = ModelResults(n).resultsPrecious.init.initSet;
    resultsPreprocessing(n) = ModelResults(n).resultsPreprocessing;
    resultsEconomic(n) = ModelResults(n).resultsEconomic;
    resultsEnvironmental(n) = ModelResults(n).resultsEnvironmental;
end

%% Plot
scatter([paramSetBase.Q],[resultsEconomic.netAnnualbeforeTax])
xlabel('Flowrate (L/s)')
ylabel('Net Annual Profit before tax (CA$)')

%% Find maximum profit simulation
idx = find([resultsEconomic.netAnnualbeforeTax] == max([resultsEconomic.netAnnualbeforeTax]));
sim_maxprofit = sim_files(idx).name;
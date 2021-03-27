%%{
clear all
%% Open Files

sim_files = dir(fullfile('FullModel', '*.mat'));
ModelResults = open(fullfile('FullModel',sim_files(1).name));
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('FullModel',sim_files(k).name));
end
%}
%% Extract Data
for n = 1:1:length(ModelResults)
    paramSetBase(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
    initSetBase(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
    paramSetPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
    initSetPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
    resultsPreprocessing(n) = ModelResults(n).ModelResults.resultsPreprocessing;
    resultsEconomic(n) = ModelResults(n).ModelResults.resultsEconomic;
    resultsEnvironmental(n) = ModelResults(n).ModelResults.resultsEnvironmental;
end

%% Annual Profit Plots
figure
subplot(2,1,1);
scatter([paramSetBase.V_app],[resultsEconomic.netAnnualbeforeTax])
xlabel('V applied (V)')
ylabel('Net Annual Profit before tax (CA$)')
subplot(2,1,2);
scatter([paramSetPrecious.V_app],[resultsEconomic.netAnnualbeforeTax])
xlabel('V applied (V)')
ylabel('Net Annual Profit before tax (CA$)')

%% Find maximum profit simulation
idx = find([resultsEconomic.netAnnualbeforeTax] == max([resultsEconomic.netAnnualbeforeTax]));
sim_maxprofit = sim_files(idx).name;
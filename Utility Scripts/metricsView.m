%%{
clear all
%% Open Files

sim_files = dir(fullfile('Simulations\FullModel3', '*.mat'));
ModelResults = open(fullfile('Simulations\FullModel3',sim_files(1).name));
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\FullModel3',sim_files(k).name));
end

%% Extract Data
for n = 1:1:length(ModelResults)
    paramSetBase(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
    initSetBase(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
    paramSetPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
    initSetPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
    resultsPreprocessing(n) = ModelResults(n).ModelResults.resultsPreprocessing;
    resultsEconomic(n) = ModelResults(n).ModelResults.resultsEconomic;
    resultsEnvironmental(n) = ModelResults(n).ModelResults.resultsEnvironmental;
    if ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod > 10 | ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod < 0 | ModelResults(n).ModelResults.resultsEnvironmental.metrics.carbonIntensity./1e6 > 0.411
        paramSetFailBase(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
        initSetFailBase(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
        paramSetFailPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
        initSetFailPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
        resultsFailPreprocessing(n) = ModelResults(n).ModelResults.resultsPreprocessing;
        resultsFailEconomic(n) = ModelResults(n).ModelResults.resultsEconomic;
        resultsFailEnvironmental(n) = ModelResults(n).ModelResults.resultsEnvironmental;
    else
        paramSetSuccessBase(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
        initSetSuccessBase(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
        paramSetSuccessPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
        initSetSuccessPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
        resultsSuccessPreprocessing(n) = ModelResults(n).ModelResults.resultsPreprocessing;
        resultsSuccessEconomic(n) = ModelResults(n).ModelResults.resultsEconomic;
        resultsSuccessEnvironmental(n) = ModelResults(n).ModelResults.resultsEnvironmental;
    end
end
%}
%% Annual Profit Plots
figure
subplot(2,1,1);
x = [initSetSuccessBase.solidPCB];
y = [resultsSuccessEconomic.metrics];
scatter(([x.m_PCB_total]./[paramSetSuccessBase.vol_bed]),[y.netAnnualafterTax])
xlabel('Loading (kg/L)')
ylabel('Net Annual Profit after tax (CA$)')
subplot(2,1,2);
x = [initSetSuccessPrecious.solidPCB];
scatter(([x.m_PCB_total]./[paramSetSuccessPrecious.vol_bed]),[y.netAnnualafterTax])
xlabel('Loading (kg/L)')
ylabel('Net Annual Profit after tax (CA$)')

%% Find maximum profit simulation
idx = find([resultsEconomic.netAnnualbeforeTax] == max([resultsEconomic.netAnnualbeforeTax]));
sim_maxprofit = sim_files(idx).name;
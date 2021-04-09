%%{
clear all
%% Open Files

sim_files = dir(fullfile('Simulations\MonteCarloStage2', '*.mat'));
ModelResults = open(fullfile('Simulations\MonteCarloStage2',sim_files(1).name));
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\MonteCarloStage2',sim_files(k).name));
end
numFail = 0;

%% Extract Data
for n = 1:1:length(ModelResults)
    paramSetBase(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
    initSetBase(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
    paramSetPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
    initSetPrecious(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
    resultsPreprocessing(n) = ModelResults(n).ModelResults.resultsPreprocessing;
    resultsPrecious(n) = ModelResults(n).ModelResults.resultsPrecious;
    resultsBase(n) = ModelResults(n).ModelResults.resultsBase;
    resultsEconomic(n) = ModelResults(n).ModelResults.resultsEconomic;
    resultsEnvironmental(n) = ModelResults(n).ModelResults.resultsEnvironmental;
    resultsBase(n) = ModelResults(n).ModelResults.resultsBase;
    if ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod > 10 | ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod < 0 | ModelResults(n).ModelResults.resultsEnvironmental.metrics.carbonIntensity./1e6 > 0.411
        paramSetBaseInfeasible(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
        initSetBaseInfeasible(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
        paramSetPreciousInfeasible(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
        initSetPreciousInfeasible(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
        resultsPreprocessingInfeasible(n) = ModelResults(n).ModelResults.resultsPreprocessing;
        resultsEconomicInfeasible(n) = ModelResults(n).ModelResults.resultsEconomic;
        resultsEnvironmentalInfeasible(n) = ModelResults(n).ModelResults.resultsEnvironmental;
        resultsBaseInfeasible(n) = ModelResults(n).ModelResults.resultsBase;
        resultsPreciousInfeasible(n) = ModelResults(n).ModelResults.resultsPrecious;
        numFail = numFail + 1;
    else
        paramSetBaseFeasible(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
        initSetBaseFeasible(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
        paramSetPreciousFeasible(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
        initSetPreciousFeasible(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
        resultsPreprocessingFeasible(n) = ModelResults(n).ModelResults.resultsPreprocessing;
        resultsEconomicFeasible(n) = ModelResults(n).ModelResults.resultsEconomic;
        resultsEnvironmentalFeasible(n) = ModelResults(n).ModelResults.resultsEnvironmental;
        resultsBaseFeasible(n) = ModelResults(n).ModelResults.resultsBase;
        resultsPreciousFeasible(n) = ModelResults(n).ModelResults.resultsPrecious;
    end
end
numSuc = length(ModelResults) - numFail;
%}
%% Annual Profit Plots
figure
subplot(2,1,1);

x = [initSetBase.solidPCB];
y = [resultsEconomic.metrics];
scatter([x.r_particles],[y.netAnnualafterTax])
xlabel('Radius (m)')
ylabel('Net Annual Profit after tax (CA$)')
%{
subplot(2,1,2);
scatter([x.r_particles],[y.netAnnualafterTax])
xlabel('Radius (m)')
ylabel('Net Annual Profit after tax (CA$)')
%}
%% Find maximum profit simulation
idx = find([resultsEconomic.netAnnualbeforeTax] == max([resultsEconomic.netAnnualbeforeTax]));
sim_maxprofit = sim_files(idx).name;
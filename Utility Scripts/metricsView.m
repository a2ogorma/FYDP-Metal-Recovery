%%{
clear all
%% Open Files

sim_files = dir(fullfile('Simulations\Sensitivity\V2.5-7.5', '*.mat'));
ModelResults = open(fullfile('Simulations\Sensitivity\V2.5-7.5',sim_files(1).name));
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\Sensitivity\V2.5-7.5',sim_files(k).name));
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
    else
        paramSetBaseFeasible(n) = ModelResults(n).ModelResults.resultsBase.init.paramSet;
        initSetBaseFeasible(n) = ModelResults(n).ModelResults.resultsBase.init.initSet;
        paramSetPreciousFeasible(n) = ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
        initSetPreciousFeasible(n) = ModelResults(n).ModelResults.resultsPrecious.init.initSet;
        resultsPreprocessingFeasible(n) = ModelResults(n).ModelResults.resultsPreprocessing;
        resultsEconomicFeasible(n) = ModelResults(n).ModelResults.resultsEconomic;
        resultsEnvironmentalFeasible(n) = ModelResults(n).ModelResults.resultsEnvironmental;
        resultsBaseFeasible(n) = ModelResults(n).ModelResults.resultsBase;
    end
end
%}
%% Annual Profit Plots
figure
subplot(2,1,1);
x = [initSetBaseFeasible.solidPCB];
y = [resultsEconomicFeasible.metrics];
%y2 = [resultsEconomicInfeasible.metrics];
scatter([paramSetBaseFeasible.V_app],[y.netAnnualafterTax],36,'b')
hold on
%scatter([resultsBaseInfeasible.numberUnits],[y2.netAnnualafterTax],36,'r')
xlabel('Applied Sensitivity\V2.5-7.5Voltage (V)')
ylabel('Net Annual Profit after tax (CA$)')
%{
subplot(2,1,2);
x = [initSetSuccessPrecious.solidPCB];
scatter(([x.m_PCB_total]./[paramSetSuccessPrecious.vol_bed]),[y.netAnnualafterTax])
xlabel('Loading (kg/L)')
ylabel('Net Annual Profit after tax (CA$)')
%}
%% Find maximum profit simulation
idx = find([resultsEconomic.netAnnualbeforeTax] == max([resultsEconomic.netAnnualbeforeTax]));
sim_maxprofit = sim_files(idx).name;
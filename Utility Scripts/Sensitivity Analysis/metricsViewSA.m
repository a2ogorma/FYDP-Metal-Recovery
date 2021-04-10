%%{
clear all
%% Open File
sim_files = dir(fullfile('Simulations\Sensitivity\p_tau_lch10-60', '*.mat'));
ModelResults = open(fullfile('Simulations\Sensitivity\p_tau_lch10-60',sim_files(1).name));
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\Sensitivity\p_tau_lch10-60',sim_files(k).name));
end

%% Extract Data
for n = 1:1:length(ModelResults.ModelResults)
    paramSetBase(n) = ModelResults.ModelResults(n).ModelResults.resultsBase.init.paramSet;
    initSetBase(n) = ModelResults.ModelResults(n).ModelResults.resultsBase.init.initSet;
    paramSetPrecious(n) = ModelResults.ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
    initSetPrecious(n) = ModelResults.ModelResults(n).ModelResults.resultsPrecious.init.initSet;
    resultsPreprocessing(n) = ModelResults.ModelResults(n).ModelResults.resultsPreprocessing;
    resultsEconomic(n) = ModelResults.ModelResults(n).ModelResults.resultsEconomic;
    resultsEnvironmental(n) = ModelResults.ModelResults(n).ModelResults.resultsEnvironmental;
    resultsBase(n) = ModelResults.ModelResults(n).ModelResults.resultsBase;
    if ModelResults.ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod > 10 | ModelResults.ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod < 0 | ModelResults.ModelResults(n).ModelResults.resultsEnvironmental.metrics.carbonIntensity > 0.411
        paramSetBaseInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsBase.init.paramSet;
        initSetBaseInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsBase.init.initSet;
        paramSetPreciousInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
        initSetPreciousInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsPrecious.init.initSet;
        resultsPreprocessingInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsPreprocessing;
        resultsEconomicInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsEconomic;
        resultsEnvironmentalInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsEnvironmental;
        resultsBaseInfeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsBase;
    else
        paramSetBaseFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsBase.init.paramSet;
        initSetBaseFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsBase.init.initSet;
        paramSetPreciousFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsPrecious.init.paramSet;
        initSetPreciousFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsPrecious.init.initSet;
        resultsPreprocessingFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsPreprocessing;
        resultsEconomicFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsEconomic;
        resultsEnvironmentalFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsEnvironmental;
        resultsBaseFeasible(n) = ModelResults.ModelResults(n).ModelResults.resultsBase;
    end
end
%}
%% Plots
idx = 1; %starting index
sol = [initSetBase.solution];
%x = [sol.Ci_Fe3_cell];
%x = [paramSetBase.V_app];
%x = [resultsPreprocessing.d_particles]*1000;
%x = [paramSetBase.tfinal]/3600;
x = [paramSetPrecious.Q];
x = x(idx:end);
met = [resultsEconomic.metrics];
met = met(idx:end);
ec = [resultsEconomic];
ec = ec(idx:end);
env = [resultsEnvironmental.metrics];
env = env(idx:end);
xaxis = 'Precious metal stage flowrate (L/s)';
figure
subplot(2,2,1);
%plot(x,[ec.totalExpenses],'b.')
plot(x, [ec.totalCapitalInvestment],'b.',x,[ec.totalExpenses],'r.',x,[met.netAnnualafterTax],'g.')
xlabel(xaxis)
ylabel('CA$')
legend('Capital Investment', 'Annual Operating Expenses', 'Net Annual Profit')
subplot(2,2,2);
plot(x,[met.paybackPeriod],'b.')
xlabel(xaxis)
ylabel('Payback period (yrs)')
subplot(2,2,3);
plot(x,[env.carbonIntensityPerMassMetal],'b.')
xlabel(xaxis)
ylabel('Carbon Intensity (t e-CO2/t metal recovered)')
subplot(2,2,4);
plot(x,[env.wasteRecoveryCorrected]*100,'b.')
xlabel(xaxis)
ylabel('Metal Recovery (%)')

%% Find maximum profit simulation
idx = find([resultsEconomic.netAnnualbeforeTax] == max([resultsEconomic.netAnnualbeforeTax]));
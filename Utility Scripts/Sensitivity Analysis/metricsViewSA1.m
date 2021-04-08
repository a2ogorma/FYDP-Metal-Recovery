%%{
clear all
%% Open File
sim_files = dir(fullfile('Simulations\Sensitivity\p_voltage3-7', '*.mat'));
ModelResults = open(fullfile('Simulations\Sensitivity\p_voltage3-7',sim_files(1).name));
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\Sensitivity\p_voltage3-7',sim_files(k).name));
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
    if ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod > 10 | ModelResults(n).ModelResults.resultsEconomic.metrics.paybackPeriod < 0 | ModelResults(n).ModelResults.resultsEnvironmental.metrics.carbonIntensity > 0.411
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
%% Plots
idx = 1; %starting index
%x = [paramSetBaseFeasible.n_units];
%sol = [initSetBaseFeasible.solution];
%x = [sol.Ci_Fe3_cell];
%x = [resultsPreprocessing.d_particles];
x = [paramSetPreciousFeasible.V_app];
%xif = [paramSetPreciousFeasible.tfinal]/3600;
x = x(idx:end);
y = [resultsEconomicFeasible.metrics];
y = y(idx:end);
ec = [resultsEconomicFeasible];
ec = ec(idx:end);
env = [resultsEnvironmentalFeasible.metrics];
env = env(idx:end);
xaxis = 'Precious metal stage voltage (V)';
figure
subplot(2,2,1);
%plot(x,[ec.totalExpenses],'b.')
plot(x,[ec.totalCapitalInvestment],'b.',x,[ec.totalExpenses],'r.')
xlabel(xaxis)
ylabel('Cost (CA$/yr)')
%ylabel('Operating Cost (CA$/yr)')
legend('Capital Investment', 'Annual Operating Expenses')
subplot(2,2,2);
plot(x,[y.netAnnualafterTax],'b.')
xlabel(xaxis)
ylabel('Net Annual Profit after tax (CA$)')
subplot(2,2,3);
plot(x,[env.carbonIntensityPerMassMetal],'b.')
xlabel(xaxis)
ylabel('Carbon Intensity (t e-CO2/t metal recovered)')
subplot(2,2,4);
plot(x,[env.wasteRecovery]*100,'b.')
xlabel(xaxis)
ylabel('Metal Recovery (%)')

%% Find maximum profit simulation
idx = find([resultsEconomic.netAnnualbeforeTax] == max([resultsEconomic.netAnnualbeforeTax]));
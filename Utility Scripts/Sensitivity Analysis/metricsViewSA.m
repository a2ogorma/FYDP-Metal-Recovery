%%{
clear all
%% Open File
sim_files = dir(fullfile('Simulations\Sensitivity\p_voltage3-7', '*.mat'));
ModelResults = open(fullfile('Simulations\Sensitivity\p_voltage3-7',sim_files(1).name));
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\Sensitivity\p_voltage3-7',sim_files(k).name));
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
idx = 2; %starting index
%x = [paramSetBaseFeasible.n_units];
%sol = [initSetBaseFeasible.solution];
%x = [sol.Ci_Fe3_cell];
x = [resultsPreprocessing.d_particles];
%x = [paramSetPreciousFeasible.n_units];
%xif = [paramSetPreciousFeasible.tfinal]/3600;
x = x(idx:end);
y = [resultsEconomicFeasible.metrics];
y = y(idx:end);
ec = [resultsEconomicFeasible];
ec = ec(idx:end);
env = [resultsEnvironmentalFeasible.metrics];
env = env(idx:end);
xaxis = 'Grinder Output Diameter (mm)';
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
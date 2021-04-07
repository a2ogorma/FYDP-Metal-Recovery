clear all
%% Open Files
sim_files = dir(fullfile('Simulations\Sensitivity\p_tau_lch10-60', '*.mat'));
ModelResults = open(fullfile('Simulations\Sensitivity\p_tau_lch10-60',sim_files(1).name));
%{
for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\Sensitivity\p_nunits2-40',sim_files(k).name));
end
%}
%% Fix precious metal flowrates and recalculate impacts
for n = 1:1:length(ModelResults.ModelResults)
    resultsPrecious = ModelResults.ModelResults(n).ModelResults.resultsPrecious;
    resultsBase = ModelResults.ModelResults(n).ModelResults.resultsBase;
    resultsPreprocessing = ModelResults.ModelResults(n).ModelResults.resultsPreprocessing;
    [resultsEnvironmental, resultsEconomic] = impactMetrics(resultsPreprocessing, resultsBase, resultsPrecious);
    ModelResults.ModelResults(n).ModelResults.resultsPrecious = resultsPrecious;
    ModelResults.ModelResults(n).ModelResults.resultsEnvironmental = resultsEnvironmental;
    ModelResults.ModelResults(n).ModelResults.resultsEconomic = resultsEconomic;
end

%% Save Files
ModelResults = ModelResults.ModelResults;
for k = 1:1:1
    save(strcat('Simulations\Sensitivity\p_tau_lch10-60\',sim_files(k).name),'ModelResults');
end
clear all
%% Open Files
sim_files = dir(fullfile('Simulations\Sensitivity\V2.5-7.5', '*.mat'));
ModelResults = open(fullfile('Simulations\Sensitivity\V2.5-7.5',sim_files(1).name));

for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\Sensitivity\V2.5-7.5',sim_files(k).name));
end

%% Fix precious metal flowrates and recalculate impacts
ModelResults = [ModelResults.ModelResults];
for n = 1:1:length(ModelResults)
    resultsPrecious = ModelResults(n).ModelResults.resultsPrecious;
    resultsBase = ModelResults(n).ModelResults.resultsBase;
    resultsPreprocessing = ModelResults(n).ModelResults.resultsPreprocessing;
    [resultsEnvironmental, resultsEconomic] = impactMetrics(resultsPreprocessing, resultsBase, resultsPrecious);
    ModelResults(n).ModelResults.resultsEnvironmental = resultsEnvironmental;
    ModelResults(n).ModelResults.resultsEconomic = resultsEconomic;
end

%% Save Files
save(strcat('Simulations\Sensitivity\V2.5-7.5\',sim_files(1).name),'ModelResults');

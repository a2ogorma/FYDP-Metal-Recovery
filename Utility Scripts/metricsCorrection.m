clear all
%% Open Files
sim_files = dir(fullfile('Simulations\FullModel3', '*.mat'));
ModelResults = open(fullfile('Simulations\FullModel3',sim_files(1).name));

for k = 2:1:length(sim_files)
    ModelResults(k) = open(fullfile('Simulations\FullModel3',sim_files(k).name));
end
%% Fix precious metal flowrates and recalculate impacts
for n = 1:1:length(ModelResults)
    resultsPrecious = ModelResults(n).ModelResults.resultsPrecious;
    resultsBase = ModelResults(n).ModelResults.resultsBase;
    resultsPreprocessing = ModelResults(n).ModelResults.resultsPreprocessing;
    [resultsEnvironmental, resultsEconomic] = impactMetrics(resultsPreprocessing, resultsBase, resultsPrecious);
    ModelResults(n).ModelResults.resultsPrecious = resultsPrecious;
    ModelResults(n).ModelResults.resultsEnvironmental = resultsEnvironmental;
    ModelResults(n).ModelResults.resultsEconomic = resultsEconomic;
end
MResults = ModelResults;
%% Save Files
for k = 1:1:length(MResults)
    ModelResults = MResults(k).ModelResults;
    save(strcat('Simulations\MonteCarloStage1\',sim_files(k).name),'ModelResults');
end
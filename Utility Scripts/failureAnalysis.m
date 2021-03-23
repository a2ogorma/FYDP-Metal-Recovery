fail = dir(fullfile('FailRuns\Precious', '*.mat'));
success = dir(fullfile('SuccessRuns\Precious','*.mat'));

result = open(fullfile('FailRuns\Precious',fail(1).name));
fail_params(1) = result.paramSetPrecious;
for k = 2:1:length(fail)
    result = open(fullfile('FailRuns\Precious',fail(k).name));
    fail_params(k) = result.paramSetPrecious;
end

result = open(fullfile('SuccessRuns\Precious',success(1).name));
success_params(1) = result.resultsPrecious.init.paramSet;

for j = 2:1:length(success)
    result = open(fullfile('SuccessRuns\Precious',success(j).name));
    success_params(j) = result.resultsPrecious.init.paramSet;
end




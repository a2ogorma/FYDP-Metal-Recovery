fail = dir(fullfile('FailPrecious', '*.mat'));
success = dir(fullfile('SuccessPrecious','*.mat'));

result = open(fullfile('FailPrecious',fail(1).name));
fail_params(1) = result.paramSetPrecious;
for k = 2:1:length(fail)
    result = open(fullfile('FailPrecious',fail(k).name));
    fail_params(k) = result.paramSetPrecious;
end

result = open(fullfile('SuccessPrecious',success(1).name));
success_params(1) = result.resultsPrecious.init.paramSet;

for j = 2:1:length(success)
    result = open(fullfile('SuccessPrecious',success(j).name));
    success_params(j) = result.resultsPrecious.init.paramSet;
end




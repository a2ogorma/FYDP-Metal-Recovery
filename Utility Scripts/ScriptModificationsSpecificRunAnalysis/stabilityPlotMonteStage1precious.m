fail = dir(fullfile('FullModel\FailedSims', '*.mat'));
success = dir(fullfile('FullModel','*.mat'));
%{
result = open(fullfile('FullModel\FailedSims',fail(1).name));
try
    fail_params(1) = result.paramSetBase;
catch
    fail_params(1) = result.resultsBase.paramSet;
end
for k = 2:1:length(fail)
    result = open(fullfile('FullModel\FailedSims',fail(k).name));
    try
        fail_params(k) = result.paramSetBase;
    catch
        fail_params(k) = result.resultsBase.paramSet;
    end
end
%}

result = open(fullfile('FullModel',success(1).name));
psuccess_params(1) = result.ModelResults.resultsPrecious.init.paramSet;
pinit_params(1) = result.ModelResults.resultsPrecious.init.initSet;
for j = 2:1:length(success)
    result = open(fullfile('FullModel',success(j).name));
    psuccess_params(j) = result.ModelResults.resultsPrecious.init.paramSet;
    pinit_params(j) = result.ModelResults.resultsPrecious.init.initSet;
end

solidPCB = [pinit_params.solidPCB];
solution = [pinit_params.solution];
r = [solidPCB.r_particles];
m = [solidPCB.m_PCB_total];
Q = [psuccess_params.Q];
Fe = [solution.Ci_Fe3_cell];
t = [psuccess_params.tfinal];
leng = [psuccess_params.length];
height = [psuccess_params.height];
n_units = [psuccess_params.n_units];
vol_bed = [psuccess_params.vol_bed];
V_app = [psuccess_params.V_app];
nbin = 20;
loading = m./vol_bed;
figure
subplot(5,2,1)
histogram(r,nbin)
title('Radius')
subplot(5,2,2)
histogram(m,nbin)
title('Mass')
subplot(5,2,3)
histogram(Q,nbin)
title('Flowrate')
subplot(5,2,4)
histogram(Fe,nbin)
title('Iron')
subplot(5,2,5)
histogram(t,nbin)
title('time')
subplot(5,2,6)
histogram(leng,nbin)
title('Cathode length')
subplot(5,2,7)
histogram(height,nbin)
title('cathode height')
subplot(5,2,8)
histogram(n_units,nbin)
title('Electrode pairs')
subplot(5,2,9)
histogram(vol_bed,nbin)
title('Bed volume')
subplot(5,2,10)
histogram(V_app,nbin)
title('Applied voltage')
sgtitle('Precious metal running values, n = 1000')
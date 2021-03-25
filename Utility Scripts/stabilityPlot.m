%run failureAnalysis before this
isSuccess = [zeros(length([fail_params.Q]),1)' ones(length([success_params.Q]),1)'];
Q = [fail_params.Q success_params.Q];
leng = [fail_params.length success_params.length];
height = [fail_params.height success_params.height];
n_units = [fail_params.n_units success_params.n_units];
vol_lch = [fail_params.vol_lch success_params.vol_lch];
tfinal = [fail_params.tfinal success_params.tfinal];

sucinit = [success_init.solidPCB];
failinit = [fail_init.solidPCB];
r_particles = [[failinit.r_particles] [sucinit.r_particles]];
m_PCB_total = [failinit.m_PCB_total sucinit.m_PCB_total];

binNo = 40;
[Qc, QbinSucRate] = barvis(Q, isSuccess, binNo);
[Lc, LbinSucRate] = barvis(leng, isSuccess, binNo);
[Hc, HbinSucRate] = barvis(height, isSuccess, binNo);
[Unitc, UnitbinSucRate] = barvis(n_units, isSuccess, binNo);
[vol_lchc, vol_lchbinSucRate] = barvis(vol_lch, isSuccess, binNo);
[tc, tbinSucRate] = barvis(tfinal, isSuccess, binNo);
[Rc, RbinSucRate] = barvis(r_particles, isSuccess, binNo);
[Mc, MbinSucRate] = barvis(m_PCB_total, isSuccess, binNo);
subplot(4,2,1)
bar(Qc,QbinSucRate)
title('Q')
subplot(4,2,2)
bar(Lc,LbinSucRate)
title('Length')
subplot(4,2,3)
bar(Hc,HbinSucRate)
title('Height')
subplot(4,2,4)
bar(Unitc,UnitbinSucRate)
title('NumberofPairs')
subplot(4,2,5)
bar(vol_lchc,vol_lchbinSucRate)
title('Leach volume')
subplot(4,2,6)
bar(tc,tbinSucRate)
title('Time')
subplot(4,2,7)
bar(Rc,RbinSucRate)
title('Radius')
subplot(4,2,8)
bar(Mc,MbinSucRate)
title('MassPCB')

function [c, binSucRate] = barvis(Q, isSuccess, binNo)
max(Q)
min(Q)
binsize = (max(Q) - min(Q))/binNo;
a = min(Q);
for j = 1:1:binNo
     b = a + binsize;
     t(j) = sum( (Q(:) >= a)&(Q(:) < b));
     Qsuc = Q.*isSuccess;
     tSuc(j) = sum( (Qsuc(:) >= a)&(Qsuc(:) < b));     
     c(j) = (a+b)/2;
     a = b;
end
binSucRate = tSuc./t;
end
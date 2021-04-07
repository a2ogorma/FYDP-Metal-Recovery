clear x
%p = anovan([resultsEconomic.totalCapitalInvestment],{[paramSetPrecious.V_app],[paramSetPrecious.Q],[paramSetPrecious.tfinal],[paramSetPrecious.length],[paramSetPrecious.height],[paramSetPrecious.n_units],[paramSetPrecious.vol_bed]},'model',1,'varnames',{'Vapp','Q','time','length','height','units','bedvol'},'continuous',[1,2,4,5,7])
    %paybackPeriod
    stats = [];
    yy = [resultsEconomic.metrics];
    y = [yy.paybackPeriod];
    stats.paybackPeriod.max = max(y);
    stats.paybackPeriod.min = min(y);
    stats.paybackPeriod.mean = mean(y);
    stats.paybackPeriod.vals = y;
    %percentOp_of_Rev
    yy = [resultsEconomic.metrics];
    y = [yy.percentOp_of_Rev];
    
    stats.percentOp_of_Rev.max = max(y);
    stats.percentOp_of_Rev.min = min(y);
    stats.percentOp_of_Rev.mean = mean(y);
    stats.percentOp_of_Rev.vals = y;
    %netAnnualafterTax
    yy = [resultsEconomic.metrics];
    y = [yy.netAnnualafterTax];
    
    stats.netAnnualafterTax.max = max(y);
    stats.netAnnualafterTax.min = min(y);
    stats.netAnnualafterTax.mean = mean(y);
    stats.netAnnualafterTax.vals = y;
    %waterIntensity
    yy = [resultsEnvironmental.metrics];
    y = [yy.waterIntensity];
    
    stats.waterIntensity.max = max(y);
    stats.waterIntensity.min = min(y);
    stats.waterIntensity.mean = mean(y);
    stats.waterIntensity.vals = y;
    %carbonIntensity
    yy = [resultsEnvironmental.metrics];
    y = [yy.carbonIntensity];
    
    stats.carbonIntensity.max = max(y);
    stats.carbonIntensity.min = min(y);
    stats.carbonIntensity.mean = mean(y);
    stats.carbonIntensity.vals = y;
    %wasteRecovery
    yy = [resultsEnvironmental.metrics];
    y = [yy.wasteRecovery];
    stats.wasteRecovery.max = max(y);
    stats.wasteRecovery.min = min(y);
    stats.wasteRecovery.mean = mean(y);
    stats.wasteRecovery.vals = y;
    %
    y = [resultsEconomic.totalCapitalInvestment];
    stats.totalCapitalInvestment.max = max(y);
    stats.totalCapitalInvestment.min = min(y);
    stats.totalCapitalInvestment.mean = mean(y);
    stats.totalCapitalInvestment.vals = y;
    %
    y = [resultsEconomic.totalRevenue];
    stats.totalRevenue.max = max(y);
    stats.totalRevenue.min = min(y);
    stats.totalRevenue.mean = mean(y);
    stats.totalRevenue.vals = y;
    
    %
    y = [resultsEconomic.totalExpenses];
    stats.totalExpenses.max = max(y);
    stats.totalExpenses.min = min(y);
    stats.totalExpenses.mean = mean(y);
    stats.totalExpenses.vals = y;
%% Fails
    %paybackPeriod
    yy = [resultsSuccessEconomic.metrics];
    y = [yy.paybackPeriod];
    statsSuc.paybackPeriod.max = max(y);
    statsSuc.paybackPeriod.min = min(y);
    statsSuc.paybackPeriod.mean = mean(y);
    statsSuc.paybackPeriod.vals = y;
    %percentOp_of_Rev
    yy = [resultsSuccessEconomic.metrics];
    y = [yy.percentOp_of_Rev];
    
    statsSuc.percentOp_of_Rev.max = max(y);
    statsSuc.percentOp_of_Rev.min = min(y);
    statsSuc.percentOp_of_Rev.mean = mean(y);
    statsSuc.percentOp_of_Rev.vals = y;
    %netAnnualafterTax
    yy = [resultsSuccessEconomic.metrics];
    y = [yy.netAnnualafterTax];
    
    statsSuc.netAnnualafterTax.max = max(y);
    statsSuc.netAnnualafterTax.min = min(y);
    statsSuc.netAnnualafterTax.mean = mean(y);
    statsSuc.netAnnualafterTax.vals = y;
    %waterIntensity
    yy = [resultsSuccessEnvironmental.metrics];
    y = [yy.waterIntensity];
    
    statsSuc.waterIntensity.max = max(y);
    statsSuc.waterIntensity.min = min(y);
    statsSuc.waterIntensity.mean = mean(y);
    statsSuc.waterIntensity.vals = y;
    %carbonIntensity
    yy = [resultsSuccessEnvironmental.metrics];
    y = [yy.carbonIntensity];
    
    statsSuc.carbonIntensity.max = max(y);
    statsSuc.carbonIntensity.min = min(y);
    statsSuc.carbonIntensity.mean = mean(y);
    statsSuc.carbonIntensity.vals = y;
    %wasteRecovery
    yy = [resultsSuccessEnvironmental.metrics];
    y = [yy.wasteRecovery];
    statsSuc.wasteRecovery.max = max(y);
    statsSuc.wasteRecovery.min = min(y);
    statsSuc.wasteRecovery.mean = mean(y);
    statsSuc.wasteRecovery.vals = y;
    %
    y = [resultsSuccessEconomic.totalCapitalInvestment];
    statsSuc.totalCapitalInvestment.max = max(y);
    statsSuc.totalCapitalInvestment.min = min(y);
    statsSuc.totalCapitalInvestment.mean = mean(y);
    statsSuc.totalCapitalInvestment.vals = y;
    %
    y = [resultsSuccessEconomic.totalRevenue];
    statsSuc.totalRevenue.max = max(y);
    statsSuc.totalRevenue.min = min(y);
    statsSuc.totalRevenue.mean = mean(y);
    statsSuc.totalRevenue.vals = y;
    
    %
    y = [resultsSuccessEconomic.totalExpenses];
    statsSuc.totalExpenses.max = max(y);
    statsSuc.totalExpenses.min = min(y);
    statsSuc.totalExpenses.mean = mean(y);
    statsSuc.totalExpenses.vals = y;
    
    yy = [resultsEconomic.metrics];
    y = [yy.valueRecovMetals];
    stats.valueRecovMetals.max = max(y);
    stats.valueRecovMetals.min = min(y);
    stats.valueRecovMetals.mean = mean(y);
    stats.valueRecovMetals.vals = y;
    
    yy = [resultsSuccessEconomic.metrics];
    y = [yy.valueRecovMetals];
    statsSuc.valueRecovMetals.max = max(y);
    statsSuc.valueRecovMetals.min = min(y);
    statsSuc.valueRecovMetals.mean = mean(y);
    statsSuc.valueRecovMetals.vals = y;
%mdl.Rsquared.Ordinary
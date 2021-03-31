clear x
%p = anovan([resultsEconomic.totalCapitalInvestment],{[paramSetPrecious.V_app],[paramSetPrecious.Q],[paramSetPrecious.tfinal],[paramSetPrecious.length],[paramSetPrecious.height],[paramSetPrecious.n_units],[paramSetPrecious.vol_bed]},'model',1,'varnames',{'Vapp','Q','time','length','height','units','bedvol'},'continuous',[1,2,4,5,7])
PsolidPCB = [initSetSuccessPrecious.solidPCB];
Psolution = [initSetSuccessPrecious.solution];
BsolidPCB = [initSetSuccessBase.solidPCB];
Bsolution = [initSetSuccessBase.solution];
x(1).name = 'Radius of initial particles (m)';
x(1).valsPrec = [PsolidPCB.r_particles];
x(1).valsBase = [BsolidPCB.r_particles];
x(2).name = 'Mass of PCB loaded (kg)';
x(2).valsPrec = [PsolidPCB.m_PCB_total];
x(2).valsBase = [BsolidPCB.m_PCB_total];
x(3).name = 'Flowrate (L/s)';
x(3).valsPrec = [paramSetSuccessPrecious.Q];
x(3).valsBase = [paramSetSuccessBase.Q];
x(4).name = 'Fe3+ conc (M)';
x(4).valsPrec = [Psolution.Ci_Fe3_cell];
x(4).valsBase = [Bsolution.Ci_Fe3_cell];
x(5).name = 'Cycle Time (h)';
x(5).valsPrec = [paramSetSuccessPrecious.tfinal]./3600;
x(5).valsBase = [paramSetSuccessBase.tfinal]./3600;
x(6).name = 'Cathode length (m)';
x(6).valsPrec = [paramSetSuccessPrecious.length];
x(6).valsBase = [paramSetSuccessBase.length];
x(7).name = 'Cathode height (m)';
x(7).valsPrec = [paramSetSuccessPrecious.height];
x(7).valsBase = [paramSetSuccessBase.height];
x(8).name = 'Number of anode-cathode pairs';
x(8).valsPrec = [paramSetSuccessPrecious.n_units];
x(8).valsBase = [paramSetSuccessBase.n_units];
x(9).name = 'Leaching vessel volume (L)';
x(9).valsPrec = [paramSetSuccessPrecious.vol_bed];
x(9).valsBase = [paramSetSuccessBase.vol_bed];
x(10).name = 'Applied voltage (V)';
x(10).valsPrec = [paramSetSuccessPrecious.V_app];
x(10).valsBase = [paramSetSuccessBase.V_app];
x(11).name = 'Mass loaded/volume (kg/L)';
x(11).valsPrec = x(2).valsPrec./x(9).valsPrec;
x(11).valsBase = x(2).valsBase./x(9).valsBase;
x(12).name = 'Total cathode area (m^2)';
x(12).valsPrec = x(6).valsPrec.*x(7).valsPrec.*x(8).valsPrec;
x(12).valsBase = x(6).valsBase.*x(7).valsBase.*x(8).valsBase;
x(13).name = 'Total cell volume (L)';
x(13).valsPrec = [paramSetSuccessPrecious.vol_cell];
x(13).valsBase = [paramSetSuccessBase.vol_cell];
x(14).name = 'Residence Time (s)';
x(14).valsPrec = x(13).valsPrec./x(3).valsPrec;
x(14).valsBase = x(13).valsBase./x(3).valsBase;
yy = [resultsSuccessEconomic.metrics];
y = [yy.netAnnualafterTax];
yval = 'Net Annual Profit after tax (CA$)';
for j = 1:1:1%size(x,2)
    xpval = x(j).valsPrec;
    xbval = x(j).valsBase;
    nameval = x(j).name;
    %netAnnualafterTax
    f = figure('DefaultAxesFontSize',15);
    subplot(2,1,1);
    scatter(xbval,y);
    xlabel(nameval, 'FontSize',17)
    ylabel(yval, 'FontSize',17)
    title('Base Metal Recovery Stage','FontSize',17)
    h = lsline;
    hslope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
    hyint = h.YData(1) - hslope*h.XData(1);
    posit = 0.25*range(xbval)+min(xbval);
    text(posit,13E5,strcat('y = ',' ',num2str(hslope,2),'*x + ',' ',num2str(hyint,2)),'FontSize',15)
    subplot(2,1,2);
    scatter(xpval,y);
    xlabel(nameval,'FontSize',17)
    ylabel(yval,'FontSize',17)
    title('Precious Metal Recovery Stage','FontSize',17)
    h = lsline;
    hslope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
    hyint = h.YData(1) - hslope*h.XData(1);
    posit = 0.25*range(xpval)+min(xpval);
    text(posit,13E5,strcat('y = ',' ',num2str(hslope,2),'*x + ',' ',num2str(hyint,2)),'FontSize',15)
    f.Position(1) = 0;
    f.Position(2) = 0;
    f.Position(3) = 1173;
    f.Position(4) = 1007;
end
%mdl.Rsquared.Ordinary
clear x
%p = anovan([resultsEconomic.totalCapitalInvestment],{[paramSetPrecious.V_app],[paramSetPrecious.Q],[paramSetPrecious.tfinal],[paramSetPrecious.length],[paramSetPrecious.height],[paramSetPrecious.n_units],[paramSetPrecious.vol_bed]},'model',1,'varnames',{'Vapp','Q','time','length','height','units','bedvol'},'continuous',[1,2,4,5,7])
PsolidPCB = [initSetSuccessPrecious.solidPCB];
Psolution = [initSetSuccessPrecious.solution];
PIsolidPCB = [initSetFailPrecious.solidPCB];
PIsolution = [initSetFailPrecious.solution];
BsolidPCB = [initSetSuccessBase.solidPCB];
Bsolution = [initSetSuccessBase.solution];
BIsolidPCB = [initSetFailBase.solidPCB];
BIsolution = [initSetFailBase.solution];
resultBase = [resultsSuccessBase];
resultPrec = [resultsSuccessPrecious];
IresultBase = [resultsFailBase];
IresultPrec = [resultsFailPrecious];
x(1).name = 'Radius of initial particles (m)';
x(1).valsFeasPrec = [PsolidPCB.r_particles];
x(1).valsFeasBase = [BsolidPCB.r_particles];
x(1).valsInPrec = [PIsolidPCB.r_particles];
x(1).valsInBase = [BIsolidPCB.r_particles];
x(1).valsAllPrec = [PsolidPCB.r_particles PIsolidPCB.r_particles];
x(1).valsAllBase = [BsolidPCB.r_particles BIsolidPCB.r_particles];
x(2).name = 'Mass of PCB loaded (kg)';
x(2).valsFeasPrec = [PsolidPCB.m_PCB_total];
x(2).valsFeasBase = [BsolidPCB.m_PCB_total];
x(2).valsInPrec = [PIsolidPCB.m_PCB_total];
x(2).valsInBase = [BIsolidPCB.m_PCB_total];
x(2).valsAllPrec = [PsolidPCB.m_PCB_total PIsolidPCB.m_PCB_total];
x(2).valsAllBase = [BsolidPCB.m_PCB_total BIsolidPCB.m_PCB_total];
x(3).name = 'Flowrate (L/s)';
x(3).valsFeasPrec = [paramSetSuccessPrecious.Q];
x(3).valsFeasBase = [paramSetSuccessBase.Q];
x(3).valsInPrec = [paramSetFailPrecious.Q];
x(3).valsInBase = [paramSetFailBase.Q];
x(3).valsAllPrec = [x(3).valsFeasPrec x(3).valsInPrec];
x(3).valsAllBase = [x(3).valsFeasBase x(3).valsInBase];
x(4).name = 'Fe3+ conc (M)';
x(4).valsFeasPrec = [Psolution.Ci_Fe3_cell];
x(4).valsFeasBase = [Bsolution.Ci_Fe3_cell];
x(4).valsInPrec = [PIsolution.Ci_Fe3_cell];
x(4).valsInBase = [BIsolution.Ci_Fe3_cell];
x(4).valsAllPrec = [x(4).valsFeasPrec x(4).valsInPrec];
x(4).valsAllBase = [x(4).valsFeasBase x(4).valsInBase];
x(5).name = 'Cycle Time (h)';
x(5).valsFeasPrec = [paramSetSuccessPrecious.tfinal]./3600;
x(5).valsFeasBase = [paramSetSuccessBase.tfinal]./3600;
x(5).valsInPrec = [paramSetFailPrecious.tfinal]./3600;
x(5).valsInBase = [paramSetFailBase.tfinal]./3600;
x(5).valsAllPrec = [x(5).valsFeasPrec x(5).valsInPrec];
x(5).valsAllBase = [x(5).valsFeasBase x(5).valsInBase];
x(6).name = 'Cathode length (m)';
x(6).valsFeasPrec = [paramSetSuccessPrecious.length];
x(6).valsFeasBase = [paramSetSuccessBase.length];
x(6).valsInPrec = [paramSetFailPrecious.length];
x(6).valsInBase = [paramSetFailBase.length];
x(6).valsAllPrec = [x(6).valsFeasPrec x(6).valsInPrec];
x(6).valsAllBase = [x(6).valsFeasBase x(6).valsInBase];
x(7).name = 'Cathode height (m)';
x(7).valsFeasPrec = [paramSetSuccessPrecious.height];
x(7).valsFeasBase = [paramSetSuccessBase.height];
x(7).valsInPrec = [paramSetFailPrecious.height];
x(7).valsInBase = [paramSetFailBase.height];
x(7).valsAllPrec = [x(7).valsFeasPrec x(7).valsInPrec];
x(7).valsAllBase = [x(7).valsFeasBase x(7).valsInBase];
x(8).name = 'Number of anode-cathode pairs';
x(8).valsFeasPrec = [paramSetSuccessPrecious.n_units];
x(8).valsFeasBase = [paramSetSuccessBase.n_units];
x(8).valsInPrec = [paramSetFailPrecious.n_units];
x(8).valsInBase = [paramSetFailBase.n_units];
x(8).valsAllPrec = [x(8).valsFeasPrec x(8).valsInPrec];
x(8).valsAllBase = [x(8).valsFeasBase x(8).valsInBase];
x(9).name = 'Leaching vessel volume (L)';
x(9).valsFeasPrec = [paramSetSuccessPrecious.vol_bed];
x(9).valsFeasBase = [paramSetSuccessBase.vol_bed];
x(9).valsInPrec = [paramSetFailPrecious.vol_bed];
x(9).valsInBase = [paramSetFailBase.vol_bed];
x(9).valsAllPrec = [x(9).valsFeasPrec x(9).valsInPrec];
x(9).valsAllBase = [x(9).valsFeasBase x(9).valsInBase];
x(10).name = 'Applied voltage (V)';
x(10).valsFeasPrec = [paramSetSuccessPrecious.V_app];
x(10).valsFeasBase = [paramSetSuccessBase.V_app];
x(10).valsInPrec = [paramSetFailPrecious.V_app];
x(10).valsInBase = [paramSetFailBase.V_app];
x(10).valsAllPrec = [x(10).valsFeasPrec x(10).valsInPrec];
x(10).valsAllBase = [x(10).valsFeasBase x(10).valsInBase];
x(11).name = 'Mass loaded/volume (kg/L)';
x(11).valsFeasPrec = x(2).valsFeasPrec./x(9).valsFeasPrec;
x(11).valsFeasBase = x(2).valsFeasBase./x(9).valsFeasBase;
x(11).valsInPrec = x(2).valsInPrec./x(9).valsInPrec;
x(11).valsInBase = x(2).valsInBase./x(9).valsInBase;
x(11).valsAllPrec = [x(11).valsFeasPrec x(11).valsInPrec];
x(11).valsAllBase = [x(11).valsFeasBase x(11).valsInBase];
x(12).name = 'Total cathode area (m^2)';
x(12).valsFeasPrec = x(6).valsFeasPrec.*x(7).valsFeasPrec.*x(8).valsFeasPrec;
x(12).valsFeasBase = x(6).valsFeasBase.*x(7).valsFeasBase.*x(8).valsFeasBase;
x(12).valsInPrec = x(6).valsInPrec.*x(7).valsInPrec.*x(8).valsInPrec;
x(12).valsInBase = x(6).valsInBase.*x(7).valsInBase.*x(8).valsInBase;
x(12).valsAllPrec = [x(12).valsFeasPrec x(12).valsInPrec];
x(12).valsAllBase = [x(12).valsFeasBase x(12).valsInBase];
x(13).name = 'Total cell volume (L)';
x(13).valsFeasPrec = [paramSetSuccessPrecious.vol_cell];
x(13).valsFeasBase = [paramSetSuccessBase.vol_cell];
x(13).valsInPrec = [paramSetFailPrecious.vol_cell];
x(13).valsInBase = [paramSetFailBase.vol_cell];
x(13).valsAllPrec = [x(13).valsFeasPrec x(13).valsInPrec];
x(13).valsAllBase = [x(13).valsFeasBase x(13).valsInBase];
x(14).name = 'Residence Time in electrowinning cell (s)';
x(14).valsFeasPrec = x(13).valsFeasPrec./x(3).valsFeasPrec;
x(14).valsFeasBase = x(13).valsFeasBase./x(3).valsFeasBase;
x(14).valsInPrec = x(13).valsInPrec./x(3).valsInPrec;
x(14).valsInBase = x(13).valsInBase./x(3).valsInBase;
x(14).valsAllPrec = [x(14).valsFeasPrec x(14).valsInPrec];
x(14).valsAllBase = [x(14).valsFeasBase x(14).valsInBase];
x(15).name = 'Residence Time in leaching unit (s)';
x(15).valsFeasPrec = x(9).valsFeasPrec./x(3).valsFeasPrec;
x(15).valsFeasBase = x(9).valsFeasBase./x(3).valsFeasBase;
x(15).valsInPrec = x(9).valsInPrec./x(3).valsInPrec;
x(15).valsInBase = x(9).valsInBase./x(3).valsInBase;
x(15).valsAllPrec = [x(15).valsFeasPrec x(15).valsInPrec];
x(15).valsAllBase = [x(15).valsFeasBase x(15).valsInBase];
x(16).name = 'Number of stage units';
x(16).valsFeasBase = [resultBase.numberUnits];
x(16).valsInBase = [IresultBase.numberUnits];
x(16).valsAllBase = [x(16).valsFeasBase x(16).valsInBase];
x(17).name = 'Leaching cell L/D ratio';
x(17).valsFeasBase = [paramSetSuccessBase.LD_bed];
x(17).valsInBase = [paramSetFailBase.LD_bed];
x(17).valsAllBase = [x(17).valsFeasBase x(17).valsInBase];
yy = [resultsSuccessEconomic.metrics];
yyI = [resultsFailEconomic.metrics];
y = [yy.netAnnualafterTax];
yI = [yyI.netAnnualafterTax];
yall = [y yI];
yval = 'Net Annual Profit after tax (CA$)';
for j = 1:1:size(x,2)
    if j == 16 | j == 17
        xbval = x(j).valsFeasBase;
        xbIval = x(j).valsInBase;
        xballval = x(j).valsAllBase;
        nameval = x(j).name;
        f = figure('DefaultAxesFontSize',15);
        scatter(xballval, yall,'MarkerEdgeColor','none');
        h = lsline;
        hslope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
        hyint = h.YData(1) - hslope*h.XData(1);
        posit = 0.75*range(xballval)+min(xballval);
        text(posit,2.7E5,strcat('y = ',' ',num2str(hslope,2),'*x + ',' ',num2str(hyint,2)),'FontSize',15)
        hold on
        scatter(xbval,y,[],'blue');
        scatter(xbIval, yI, [], 'red');
        legend('Legend','Trendline','Feasible runs','Infeasible runs')
        hold off
        xlabel(nameval, 'FontSize',17)
        ylabel(yval, 'FontSize',17)
        title('Base Metal Recovery Stage','FontSize',17)
        f.Position(1) = 0;
        f.Position(2) = 0;
        f.Position(3) = 1173;
        f.Position(4) = 607;
    else
        xpval = x(j).valsFeasPrec;
        xpIval = x(j).valsInPrec;
        xbval = x(j).valsFeasBase;
        xbIval = x(j).valsInBase;
        xpallval = x(j).valsAllPrec;
        xballval = x(j).valsAllBase;
        nameval = x(j).name;
        %netAnnualafterTax
        f = figure('DefaultAxesFontSize',15);
        subplot(2,1,1);
        scatter(xballval, yall,'MarkerEdgeColor','none');
        h = lsline;
        hslope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
        hyint = h.YData(1) - hslope*h.XData(1);
        posit = 0.75*range(xballval)+min(xballval);
        text(posit,2.7E5,strcat('y = ',' ',num2str(hslope,2),'*x + ',' ',num2str(hyint,2)),'FontSize',15)
        hold on
        scatter(xbval,y, [], 'blue');
        scatter(xbIval, yI, [], 'red');
        legend('Legend','Trendline','Feasible runs','Infeasible runs')
        hold off
        xlabel(nameval, 'FontSize',17)
        ylabel(yval, 'FontSize',17)
        title('Base Metal Recovery Stage','FontSize',17)
        subplot(2,1,2);
        scatter(xpallval, yall,'MarkerEdgeColor','none');
        h = lsline;
        hslope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
        hyint = h.YData(1) - hslope*h.XData(1);
        posit = 0.75*range(xpallval)+min(xpallval);
        text(posit,2.7E5,strcat('y = ',' ',num2str(hslope,2),'*x + ',' ',num2str(hyint,2)),'FontSize',15)
        hold on
        scatter(xpval,y, [], 'blue');    
        scatter(xpIval, yI, [], 'red');
        legend('Legend','Trendline','Feasible runs','Infeasible runs')
        hold off
        xlabel(nameval,'FontSize',17)
        ylabel(yval,'FontSize',17)
        title('Precious Metal Recovery Stage','FontSize',17)
        f.Position(1) = 0;
        f.Position(2) = 0;
        f.Position(3) = 1173;
        f.Position(4) = 1007;
    end
end
%mdl.Rsquared.Ordinary
%}
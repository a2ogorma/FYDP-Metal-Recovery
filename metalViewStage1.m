%Plots for stage one extraction/recovery of base metals
%Select metal here
metal = 1;
metal_names = {'Copper';'Tin';'Iron';'Silver';'Gold';'Palladium'};
ion_names = {'Cu2+', 'Sn2+', 'Fe2+', 'Fe3+', 'Ag+', 'Au3+', 'Pd2+'};
propertiesMetals
results = resultsBase;
init = initSetBase;
param = paramSetBase;
vol_cell = param.vol_cell;
vol_lch = param.vol_lch;
t = results.t;
Cm = results.Cm;
Erev_cat = results.electrowinning.Erev_cat;
Erev_an = results.electrowinning.Erev_an;
Erev_lch = results.leaching.Erev_lch;
E_corr = results.leaching.E_corr;
E_cat = results.electrowinning.E_cat;
E_an = results.electrowinning.E_an;
I_cell = results.electrowinning.I_cell;
I_corr = results.leaching.I_corr;
m_plated = results.electrowinning.m_plated;
m_plated = m_plated - m_plated(1,:);
m_PCB = results.PCB.massRem(:,2:7);
m_PCB_i = m_PCB(1,:);
m_lch = m_PCB_i - m_PCB;
pct_lch = 100*m_lch./m_PCB_i;
m_sol_i = ([Cm(1,1:2) Cm(1,3)+Cm(1,4) Cm(1,5:7)]*...
    vol_cell+[Cm(1,21:22) Cm(1,23)+Cm(1,24) Cm(1,25:27)]*vol_lch).*mw(2:7)/1000; %mass of metals, kg in solution at t=0
pct_rec = subplus(100*(m_plated-m_sol_i)./m_PCB_i);
f = figure;
set(f, 'DefaultLegendLocation', 'southwest');
sgtitle(metal_names(metal))
set(gcf, 'Position',  [40, 40, 1500, 700])

%find where solver failed and eliminate those time steps
te = results.te;
ie = results.ie;
for j = 1:1:numel(te)
    t_index = find(te(j));
end

t = t/3600; %Convert to hrs
i = 500; %initial index
tf = size(t);
if metal == 3 %Iron -- Two ions in this case
    Fe2 = 3;
    Fe3 = 4;
    rxn1 = 3;
    rxn2 = 4;
    %Concentration Plot Fe2+
    subplot(2,4,1);
    plot(t(i:tf),Cm(i:tf,Fe2),t(i:tf),Cm(i:tf,Fe2+10),t(i:tf),Cm(i:tf,Fe2+20));
    title('Fe2+ Concentration');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (hrs.)');
    ylabel('Concentration (M)');
    %Concentration Plot Fe2+
    subplot(2,4,2);
    plot(t(i:tf),Cm(i:tf,Fe3),t(i:tf),Cm(i:tf,Fe3+10),t(i:tf),Cm(i:tf,Fe3+20));
    title('Fe3+ Concentration');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (hrs.)');
    ylabel('Concentration (M)');
    %Current for Fe3+ + e- <--> Fe2+
    subplot(2,4,3);
    plot(t(i:tf),I_cell(i:tf,rxn1),t(i:tf),I_corr(i:tf,rxn1));
    title('Current: Fe3+/Fe2+')
    legend('Electrowinning','Leaching')
    xlabel('Time (hrs.)')
    ylabel('Current (A)')
    %Current for Fe2+ + e- <--> Fe(s)
    subplot(2,4,4);
    plot(t(i:tf),I_cell(i:tf,rxn2),t(i:tf),I_corr(i:tf,rxn2));
    title('Current: Fe2+/Fe(s)')
    legend('Electrowinning','Leaching')
    xlabel('Time (hrs.)')
    ylabel('Current (A)')
    %Mass
    subplot(2,4,5);
    plot(t(i:tf),pct_lch(i:tf,metal),t(i:tf),pct_rec(i:tf,metal));
    title(['Initial mass: ' num2str(m_PCB_i(metal)) ' kg']);
    legend('Metal leached','Metal recovered');
    xlabel('Time (hrs.)');
    ylabel('% mass');
    %Nernst
    subplot(2,4,6);
    plot(t(i:tf),Erev_cat(i:tf,rxn1),t(i:tf),Erev_cat(i:tf,rxn2),t(i:tf),E_cat(i:tf));
    title('Potentials');
    legend('Fe3+/Fe2+ Erev in Catholyte', 'Fe2+/Fe(s) Erev in Catholyte', 'Ecat');
    subplot(2,4,7);
    plot(t(i:tf),Erev_an(i:tf,rxn1),t(i:tf),Erev_an(i:tf,rxn2),t(i:tf),E_an(i:tf));
    title('Potentials');
    legend('Fe3+/Fe2+ Erev in Anolyte', 'Fe2+/Fe(s) Erev in Anolyte', 'Ean');
    subplot(2,4,8);
    plot(t(i:tf),Erev_lch(i:tf,rxn1),t(i:tf),Erev_lch(i:tf,rxn2),t(i:tf),E_corr(i:tf));
    title('Potentials');
    legend('Fe3+/Fe2+ Erev in Leach unit', 'Fe2+/Fe(s) Erev in Leach unit', 'Ecorr');
    
elseif metal > 3 %Precious metals
    cation = metal+1;
    rxn = 2*cation-5;
    %Concentration Plot
    subplot(2,3,1);
    plot(t(i:tf),Cm(i:tf,cation),t(i:tf),Cm(i:tf,cation+10),t(i:tf),Cm(i:tf,cation+20));
    title('Concentrations');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (hrs.)');
    ylabel('Concentration (M)');
    %Mass
    subplot(2,3,3);
    plot(t(i:tf),pct_lch(i:tf,metal),t(i:tf),pct_rec(i:tf,metal));
    title(['Initial mass: ' num2str(m_PCB_i(metal)) ' kg']);
    legend('Metal leached','Metal recovered');
    xlabel('Time (hrs.)');
    ylabel('% mass');
    if metal == 4 %Ag
        rxnCl = rxn+1;
        %Currents
        subplot(2,3,2);
        plot(t(i:tf),I_cell(i:tf,rxn),t(i:tf),I_corr(i:tf,rxn)+I_corr(i:tf,rxnCl));
        title('Currents')
        legend('Electrowinning','Leaching')
        xlabel('Time (hrs.)')
        ylabel('Current (A)')
        %Nernst
        subplot(2,3,4);
        plot(t(i:tf),Erev_cat(i:tf,rxn),t(i:tf),Erev_cat(i:tf,rxnCl),t(i:tf),E_cat(i:tf));
        title('Catholyte Potentials');
        legend('Erev, Ag+/Ag(s)', 'Erev, AgCl/Ag(s)/Cl-', 'Ecat')
        subplot(2,3,5);
        plot(t(i:tf),Erev_an(i:tf,rxn),t(i:tf),Erev_an(i:tf,rxnCl),t(i:tf),E_an(i:tf));
        title('Anolyte Potentials');
        legend('Erev, Ag+/Ag(s)', 'Erev, AgCl/Ag(s)/Cl-', 'Ean');
        subplot(2,3,6);
        plot(t(i:tf),Erev_lch(i:tf,rxn),t(i:tf),Erev_lch(i:tf,rxnCl),t(i:tf),E_corr(i:tf));
        title('Leaching Unit Potentials');
        legend('Erev, Ag+/Ag(s)', 'Erev, AgCl/Ag(s)/Cl-', 'Ecorr')
    elseif metal == 5 %Au
        rxnCl = rxn+1;
        %Currents
        subplot(2,3,2);
        plot(t(i:tf),I_cell(i:tf,rxn),t(i:tf),I_corr(i:tf,rxn)+I_corr(i:tf,rxnCl));
        title('Currents')
        legend('Electrowinning','Leaching')
        xlabel('Time (hrs.)')
        ylabel('Current (A)')
        %Nernst
        subplot(2,3,4);
        plot(t(i:tf),Erev_cat(i:tf,rxn),t(i:tf),Erev_cat(i:tf,rxnCl),t(i:tf),E_cat(i:tf));
        title('Catholyte Potentials');
        legend('Erev, Au3+/Au(s)', 'Erev, AuCl4-/Au(s)/Cl-', 'Ecat')
        subplot(2,3,5);
        plot(t(i:tf),Erev_an(i:tf,rxn),t(i:tf),Erev_an(i:tf,rxnCl),t(i:tf),E_an(i:tf));
        title('Anolyte Potentials');
        legend('Erev, Au3+/Au(s)', 'Erev, AuCl4-/Au(s)/Cl-', 'Ean');
        subplot(2,3,6);
        plot(t(i:tf),Erev_lch(i:tf,rxn),t(i:tf),Erev_lch(i:tf,rxnCl),t(i:tf),E_corr(i:tf));
        title('Leaching Unit Potentials');
        legend('Erev, Au3+/Au(s)', 'Erev, AuCl4-/Au(s)/Cl-', 'Ecorr')
    else %Pd
        %Currents
        subplot(2,3,2);
        plot(t(i:tf),I_cell(i:tf,rxn),t(i:tf),I_corr(i:tf,rxn));
        title('Currents')
        legend('Electrowinning','Leaching')
        xlabel('Time (hrs.)')
        ylabel('Current (A)')
        %Nernst
        subplot(2,3,4);
        plot(t(i:tf),Erev_cat(i:tf,rxn),t(i:tf),E_cat(i:tf));
        title('Catholyte Potentials');
        legend('Erev, Pd2+/Pd(s)', 'Ecat')
        subplot(2,3,5);
        plot(t(i:tf),Erev_an(i:tf,rxn),t(i:tf),E_an(i:tf));
        title('Anolyte Potentials');
        legend('Erev, Pd2+/Pd(s)', 'Ean');
        subplot(2,3,6);
        plot(t(i:tf),Erev_lch(i:tf,rxn),t(i:tf),E_corr(i:tf));
        title('Leaching Unit Potentials');
        legend('Erev, Pd2+/Pd(s)', 'Ecorr')
    end
else
    cation = metal;
    rxn = metal;
    
    %Concentration Plot
    subplot(2,3,1);
    plot(t(i:tf),Cm(i:tf,cation),t(i:tf),Cm(i:tf,cation+10),t(i:tf),Cm(i:tf,cation+20));
    title('Concentrations');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (hrs.)');
    ylabel('Concentration (M)');
    %Currents
    subplot(2,3,2);
    plot(t(i:tf),I_cell(i:tf,rxn),t(i:tf),I_corr(i:tf,rxn));
    title('Currents')
    legend('Electrowinning','Leaching')
    xlabel('Time (hrs.)')
    ylabel('Current (A)')
    %Mass
    subplot(2,3,3);
    plot(t(i:tf),pct_lch(i:tf,metal),t(i:tf),pct_rec(i:tf,metal));
    title(['Initial mass: ' num2str(m_PCB_i(metal)) ' kg']);
    legend('Metal leached','Metal recovered');
    xlabel('Time (hrs.)');
    ylabel('% mass');
    %Nernst
    subplot(2,3,4);
    plot(t(i:tf),Erev_cat(i:tf,rxn),t(i:tf),E_cat(i:tf));
    title('Potentials');
    legend('Rev. Potential in Catholyte', 'Cathode Potential');
    subplot(2,3,5);
    plot(t(i:tf),Erev_an(i:tf,rxn),t(i:tf),E_an(i:tf));
    title('Potentials');
    legend('Rev. Potential in Anolyte', 'Anode Potential');
    subplot(2,3,6);
    plot(t(i:tf),Erev_lch(i:tf,rxn),t(i:tf),E_corr(i:tf));
    title('Potentials');
    legend('Rev. Potential in Leaching Unit', 'Corrosion Potential');
end
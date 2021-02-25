%Select metal here
metal = 1;
metal_names = {'Copper';'Tin';'Iron';'Silver';'Gold';'Palladium'};
ion_names = {'Cu2+', 'Sn2+', 'Fe2+', 'Fe3+', 'Ag+', 'Au3+', 'Pd2+'};
propertiesMetals
results = resultsBase;
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
m_PCB = results.PCB.massRem(:,2:7);
m_PCB_i = m_PCB(1,:);
pct_lch = 100*m_PCB./m_PCB_i;
pct_rec = 100*m_plated./m_PCB_i;
f = figure;
set(f, 'DefaultLegendLocation', 'southwest');
sgtitle(metal_names(metal))
set(gcf, 'Position',  [40, 40, 1500, 700])

if metal == 3 %Iron -- Two ions in this case
    Fe2 = 3;
    Fe3 = 4;
    rxn1 = 3;
    rxn2 = 4;
    %Concentration Plot Fe2+
    subplot(2,4,1);
    plot(t,Cm(:,Fe2),t,Cm(:,Fe2+10),t,Cm(:,Fe2+20));
    title('Fe2+ Concentration');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (s)');
    ylabel('Concentration (M)');
    %Concentration Plot Fe2+
    subplot(2,4,2);
    plot(t,Cm(:,Fe3),t,Cm(:,Fe3+10),t,Cm(:,Fe3+20));
    title('Fe3+ Concentration');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (s)');
    ylabel('Concentration (M)');
    %Current for Fe3+ + e- <--> Fe2+
    subplot(2,4,3);
    plot(t,I_cell(:,rxn1),t,I_corr(:,rxn1));
    title('Current: Fe3+/Fe2+')
    legend('Electrowinning','Leaching')
    xlabel('Time (s)')
    ylabel('Current (A)')
    %Current for Fe2+ + e- <--> Fe(s)
    subplot(2,4,4);
    plot(t,I_cell(:,rxn2),t,I_corr(:,rxn2));
    title('Current: Fe2+/Fe(s)')
    legend('Electrowinning','Leaching')
    xlabel('Time (s)')
    ylabel('Current (A)')
    %Mass
    subplot(2,4,5);
    plot(t,pct_lch(:,metal),t,pct_rec(:,metal));
    title(['Initial mass: ' num2str(m_PCB_i(metal)) ' kg']);
    legend('Metal leached','Metal recovered');
    xlabel('Time (s)');
    ylabel('% mass');
    %Nernst
    subplot(2,4,6);
    plot(t,Erev_cat(:,rxn1),t,Erev_cat(:,rxn2),t,E_cat);
    title('Potentials');
    legend('Fe3+/Fe2+ Erev in Catholyte', 'Fe2+/Fe(s) Erev in Catholyte', 'Ecat');
    subplot(2,4,7);
    plot(t,Erev_an(:,rxn1),t,Erev_an(:,rxn2),t,E_an);
    title('Potentials');
    legend('Fe3+/Fe2+ Erev in Anolyte', 'Fe2+/Fe(s) Erev in Anolyte', 'Ean');
    subplot(2,4,8);
    plot(t,Erev_lch(:,rxn1),t,Erev_lch(:,rxn2),t,E_corr);
    title('Potentials');
    legend('Fe3+/Fe2+ Erev in Leach unit', 'Fe2+/Fe(s) Erev in Leach unit', 'Ecorr');
    
elseif metal > 3 %Precious Metals
    cation = metal+1;
    rxn = 2*cation-5;
    %Concentration Plot
    subplot(2,3,1);
    plot(t,Cm(:,cation),t,Cm(:,cation+10),t,Cm(:,cation+20));
    title('Concentrations');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (s)');
    ylabel('Concentration (M)');
    %Currents
    subplot(2,3,2);
    plot(t,I_cell(:,rxn),t,I_corr(:,rxn));
    title('Currents')
    legend('Electrowinning','Leaching')
    xlabel('Time (s)')
    ylabel('Current (A)')
    %Mass
    subplot(2,3,3);
    plot(t,pct_lch(:,metal),t,pct_rec(:,metal));
    title(['Initial mass: ' num2str(m_PCB_i(metal)) ' kg']);
    legend('Metal leached','Metal recovered');
    xlabel('Time (s)');
    ylabel('% mass');
    %Nernst
    subplot(2,3,4);
    plot(t,Erev_cat(:,rxn),t,E_cat);
    title('Potentials');
    legend('Rev. Potential in Catholyte', 'Cathode Potential');
    subplot(2,3,5);
    plot(t,Erev_an(:,rxn),t,E_an);
    title('Potentials');
    legend('Rev. Potential in Anolyte', 'Anode Potential');
    subplot(2,3,6);
    plot(t,Erev_lch(:,rxn),t,E_corr);
    title('Potentials');
    legend('Rev. Potential in Leaching Unit', 'Corrosion Potential');
else
    cation = metal;
    rxn = metal;
    
    %Concentration Plot
    subplot(2,3,1);
    plot(t,Cm(:,cation),t,Cm(:,cation+10),t,Cm(:,cation+20));
    title('Concentrations');
    legend('Catholyte','Anolyte', 'Leaching');
    xlabel('Time (s)');
    ylabel('Concentration (M)');
    %Currents
    subplot(2,3,2);
    plot(t,I_cell(:,rxn),t,I_corr(:,rxn));
    title('Currents')
    legend('Electrowinning','Leaching')
    xlabel('Time (s)')
    ylabel('Current (A)')
    %Mass
    subplot(2,3,3);
    plot(t,pct_lch(:,metal),t,pct_rec(:,metal));
    title(['Initial mass: ' num2str(m_PCB_i(metal)) ' kg']);
    legend('Metal leached','Metal recovered');
    xlabel('Time (s)');
    ylabel('% mass');
    %Nernst
    subplot(2,3,4);
    plot(t,Erev_cat(:,rxn),t,E_cat);
    title('Potentials');
    legend('Rev. Potential in Catholyte', 'Cathode Potential');
    subplot(2,3,5);
    plot(t,Erev_an(:,rxn),t,E_an);
    title('Potentials');
    legend('Rev. Potential in Anolyte', 'Anode Potential');
    subplot(2,3,6);
    plot(t,Erev_lch(:,rxn),t,E_corr);
    title('Potentials');
    legend('Rev. Potential in Leaching Unit', 'Corrosion Potential');
    
end
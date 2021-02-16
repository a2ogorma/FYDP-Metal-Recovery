metal = 1
names = {'Copper'; 'Tin';'Aluminum';'Lead';'IronII';'IronIII';'Gold';'Silver';'Palladium';'Hydrogen';'Chloride'};
figure
sgtitle(names(metal))
    
subplot(2,2,1)
plot(t,Cm(:,metal),t,Cm(:,metal+11))
title('Concentrations')
legend('Electrowinning','Leaching')
xlabel('Time (s)')
ylabel('Concentrations (M)')

subplot(2,2,2)
plot(t,Erev_cell(:,metal),t,Erev_lch(:,metal),t,E_corr)
title('Nernst Potentials')
legend('Electrowinning','Leaching','Overall corrosion')
xlabel('Time (s)')
ylabel('Potential (V)')

subplot(2,2,3)
plot(t,I_cell(:,metal),t,I_corr(:,metal))
title('Currents')
legend('Electrowinning','Leaching')
xlabel('Time (s)')
ylabel('Current (A)')

%Solid mass order: Inert Cu Sn Al Pb Fe Ag Au Pd
massRem = m_PCB(:,metal+1);
massStart = m_PCB(1,metal+1);
pctRem = massRem./massStart;
massRec(1) = 0;
for j = 2:1:length(I_cell(:,metal))
    massRec(j) = -mw(metal+1)*trapz(t(1:j), I_cell(1:j,metal))/F/z(metal)/1000;
end
pctRec = -massRec./massStart;
subplot(2,2,4)
plot(t,pctRem,t,pctRec)
title(strcat('Remaining Mass (',num2str(massRem(end),2),' kg) / Mass Recovered (',num2str(massRec(end),2),' kg)'))
legend('Remaining solids','Recovered solids')
xlabel('Time (s)')
ylabel('Percent (%)')
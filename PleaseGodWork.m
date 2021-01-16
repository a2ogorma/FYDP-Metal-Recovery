clear all
%inlet concs for now
CAg0(1) = 0.10;
CAu0(1) = 0.05;
CPd0(1) = 0.0001;
CH0(1) = 1e-10;
CFe2 = 0.1;%temp
CFe3 = 0.1;%temp

%conductivity values, check on this (S m2 /mol)
gammaAg = 61.90e-4;
gammaAu = 4.1e7;
gammaPd = 0.1;
gammaH = 349.81e-4;
%gammaFe2
%gammeFe3

%divergence from ideality values (check on this shit too but low
%priority)... something about fugacity comes to mind???
gamFe2 = 1;
gamFe3 = 1;
gamAg = 1;%assume all ideal unless proven otherwise
gamAu = 1;
gamPd = 1;
gamH = 1;

%Universal constants (keep em handy)
F = 96500; %C/mol
Rgas = 8.3145; %J/molK
mmAg = 196.96657; %g/mol
mmAu = 107.8682; %g/mol
mmPd = 106.42; %g/mol

%System variables
Fin = 10; %L/s
Fout = Fin; %kept same for now, may want to make into different values to account for accumulation and controls
T = 298.15; %K, we can play around with this but if we want to vary this kinetically then I oop
V = 100; %L
cursivel = 10; %m, characteristic distance
Acat = 1; %m2, area
Aan = 10; %m2, area
Vapp = 10; %V
Rhardware = 0; %set to 0 for now before i do something with it

%setting up initial conc in cell
nAg0(1) = CAg0(1)*Fin; %mol/s
nAu0(1) = CAu0(1)*Fin; %mol/s
nPd0(1) = CPd0(1)*Fin; %mol/s
nH0(1) = CH0(1)*Fin; %mol/s
nAg(1) = CAg0(1)*V; % moles of silver remaining
nAu(1) = CAu0(1)*V; % moles
nPd(1) = CPd0(1)*V; % moles
nH(1) = CH0(1)*V; % moles

aO2 = 0.21;%atm, its always this cause atmosphere,but maybe we wanna pressurize
aH2 = 0.21;%atm - not right so FIX THIS LATER

%calculated vars (will be iterated)
K(1) = gammaAg*nAg(1)+gammaAu*nAu(1)+gammaPd*2*nPd(1)+gammaH*nH(1);
R(1) = cursivel/(Acat*K(1));
ErevAn(1) = -1.23 + (Rgas*T/(4*F))*log(aO2*(gamH*nH(1)/V)^4);
ErevFe(1) = -0.77 + (Rgas*T/F)*log(gamFe3*CFe3/(gamFe2*CFe2));
ErevAg(1) = 0.7989 + (Rgas*T/(F))*log(1/(gamAg*nAg(1)/V));
ErevAu(1) = 1.692 + (Rgas*T/(F))*log(1/(gamAg*nAu(1)/V));
ErevPd(1) = 0.951 + (Rgas*T/(2*F))*log(1/(((gamAg*nPd(1))/V)^2));
%ErevH(1) = (Rgas*T/F)*log(

%assume well mixed, or bulk and surface concentrations are the same
alphaAn = 0.5;
alphaFe = 0.5;
alphaAg = 0.5; %base assumption as most systems are close to this number
alphaAu = 0.5; %base assumption as most systems are close to this number
alphaPd = 0.5; %base assumption as most systems are close to this number
iAn0 = 1E-12; %A/m2
iFe0 = 10^-6.8;%from course notes using palladium electrode, E-4 used to convert from A/cm2 to A/m2
iAg0 = 7.3e-9; %A/m2, from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1007/s10800-007-9434-x - verify this link applies
iAu0 = 2E-8; %A/m2, value of dissolution including Na2S in dissolution reaction from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1134/S1023193506040021, converting A/cm2 to A/m2 
iPd0 = 4e-9; %using a mix of a few from that reference
iH0 = 1E-10; %need to identify suitable electrode to help limit hydrogen evolution. From volcano plot in class, likely copper, gold or silver electrodes
etaAn0 = 0.8;%initial guesses for overpot
etaCat0 = -0.8;%initial guesses for overpot

CurrentFunc = @(eta) Aan*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T))) + iFe0*(exp(alphaFe*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))-exp(-(1-alphaFe)*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))))+Acat*( iAg0*(exp(alphaAg*F*eta(2)/(Rgas*T))-exp(-(1-alphaAg)*F*eta(2)/(Rgas*T))) + iAu0*(exp(alphaAu*F*(ErevAg(1)+eta(2)-ErevAu(1))/(Rgas*T))-exp(-(1-alphaAu)*F*(ErevAg(1)+eta(2)-ErevAu(1))/(Rgas*T))) + iPd0*(exp(alphaPd*2*F*(ErevAg(1)+eta(2)-ErevPd(1))/(Rgas*T))-exp(-(1-alphaPd)*2*F*(ErevAg(1)+eta(2)-ErevPd(1))/(Rgas*T))));
VappFunc = @(eta) abs(ErevAg(1)-ErevAn(1))+Aan*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T)))+ iFe0*(exp(alphaFe*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))-exp(-(1-alphaFe)*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))))*(R(1)+Rhardware)+abs(eta(1))+abs(eta(2))-Vapp;
overpotentialFunc = @(eta) [CurrentFunc(eta);VappFunc(eta)];  
etaGuess0 = [etaAn0 etaCat0];
options = optimoptions('fsolve','MaxFunctionEvaluations',1000,'MaxIterations',1000);
eta = fsolve(overpotentialFunc,etaGuess0,options);%,[0,-Vapp],[Vapp,0]);
etaAn = eta(1);
etaFe = -ErevAn(1)+eta(1)+ErevFe;
etaAg = eta(2);
etaAu = ErevAg(1)+eta(2)-ErevAu(1);
etaPd = ErevAg(1)+eta(2)-ErevPd(1);
disp([etaAn,etaFe,etaAg,etaAu,etaPd])
CurrentFunc = Aan*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T))) + iFe0*(exp(alphaFe*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))-exp(-(1-alphaFe)*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))))+Acat*( iAg0*(exp(alphaAg*F*eta(2)/(Rgas*T))-exp(-(1-alphaAg)*F*eta(2)/(Rgas*T))) + iAu0*(exp(alphaAu*F*(ErevAg(1)+eta(2)-ErevAu(1))/(Rgas*T))-exp(-(1-alphaAu)*F*(ErevAg(1)+eta(2)-ErevAu(1))/(Rgas*T))) + iPd0*(exp(alphaPd*2*F*(ErevAg(1)+eta(2)-ErevPd(1))/(Rgas*T))-exp(-(1-alphaPd)*2*F*(ErevAg(1)+eta(2)-ErevPd(1))/(Rgas*T))))
VappFunc = abs(ErevAg(1)-ErevAn(1))+Aan*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T)))+ iFe0*(exp(alphaFe*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))-exp(-(1-alphaFe)*F*(ErevAn(1)+eta(1)-ErevFe(1))/(Rgas*T))))*(R(1)+Rhardware)+abs(eta(1))+abs(eta(2))-Vapp 
ErevAg(1)+eta(2)
ErevAn(1)+eta(1)
clear all
tic
%inlet concs for now
CAgS2O320(1) = 0.10;
CAuS2O320(1) = 0.05;
CPdS2O340(1) = 0.0001;
CH0(1) = 1e-10;
CFe20(1) = 0.005;
CFe30(1) = 0.001;
CS2O30(1) = 0.05;
COH0(1) = (1e-14)/CH0(1);

%conductivity values, check on this (S m2 /mol)
gammaAg = 61.90e-4;
gammaAu = 4.1e7;
gammaPd = 0.1;
gammaH = 349.81e-4;
gammaFe2 = 0;%easy update
gammaFe3 = 0;%easy update
gammaOH = 0;%easy update
%note many conductivity values are missing due to lack of immediately
%available data. This is sorta ok because as it is this makes a more
%conservative (aka worse/higher) resistance estimate, and it is not
%expected that the complexes would contribute as much to the system.
%Ideally, thiosulfate conductivity would be found which would allow for
%estimation of it if they werent in their complexes.

%divergence from ideality values (check on this shit too but low
%priority)... something about fugacity comes to mind???
gamAg = 1;%assume all ideal unless proven otherwise
gamAu = 1;
gamPd = 1;
gamH = 1;
gamFe2 = 1;
gamFe3 = 1;
gamOH = 1;
gamS2O3 = 1;

%Universal constants (keep em handy)
F = 96485; %C/mol
Rgas = 8.3145; %J/molK
mmAg = 196.96657; %g/mol
mmAu = 107.8682; %g/mol
mmPd = 106.42; %g/mol

%System parameters
Fin = 0; %L/s
Fout = Fin; %kept same for now, may want to make into different values to account for accumulation and controls
T = 298.15; %K, we can play around with this but if we want to vary this kinetically then I oop
V = 100; %L
cursivel = 10; %m, characteristic distance
Acat = 10; %m2, area
Aan = 10; %m2, area
Vapp = 1; %V
Rhardware = 0; %set to 0 for now before i do something with it
tfinal = 72;

%setting up initial conc in cell
nAgS2O32(1) = CAgS2O320(1)*V; % moles of silver remaining
nAuS2O32(1) = CAuS2O320(1)*V; % moles
nPdS2O34(1) = CPdS2O340(1)*V; % moles
nS2O3(1) = CS2O30(1)*V;
nFe2(1) = CFe20(1)*V;
nFe3(1) = CFe30(1)*V;
nH(1) = CH0(1)*V; % moles
nOH(1) = COH0(1)*V;
aO2 = 0.21;%atm, its always this cause atmosphere,but maybe we wanna pressurize
aH2 = 0.0001;%dunnno what this should be

%assume well mixed, or bulk and surface concentrations are the same
alphaAn = 0.5;
alphaFe = 0.5;
alphaAg = 0.5; %base assumption as most systems are close to this number
alphaAu = 0.5; %base assumption as most systems are close to this number
alphaPd = 0.5; %base assumption as most systems are close to this number
alphaH = 0.5;

iAn0 = 1E-12; %A/m2
iFe0 = 10^-6.8;%from course notes using palladium electrode, E-4 used to convert from A/cm2 to A/m2
iAg0 = 7.3e-9; %A/m2, from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1007/s10800-007-9434-x - verify this link applies
iAu0 = 2E-8; %A/m2, value of dissolution including Na2S in dissolution reaction from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1134/S1023193506040021, converting A/cm2 to A/m2 
iPd0 = 4e-9; %using a mix of a few from that reference
iH0 = 1E-10; %need to identify suitable electrode to help limit hydrogen evolution. From volcano plot in class, likely copper, gold or silver electrodes

etaAn(1) = 0.8;%initial guesses for overpot
etaAg(1) = -0.8;%initial guesses for overpot

%electrochemical equations and relations
dnAgdt = @(CAg0,nAg,t,iAg) CAg0*Fin-nAg*Fout/V+(1/F)*iAg*Acat;
dnAudt = @(CAu0,nAu,t,iAu) CAu0*Fin-nAu*Fout/V+(1/F)*iAu*Acat;
dnPddt = @(CPd0,nPd,t,iPd) CPd0*Fin-nPd*Fout/V+(1/(2*F))*iPd*Acat;
dnFe2dt = @(CFe20,nFe2,t,iFe2) CFe20*Fin-nFe2*Fout/V+(1/F)*iFe2*Aan;
dnFe3dt = @(CFe30,nFe3,t,iFe2) CFe30*Fin-nFe3*Fout/V-(1/F)*iFe2*Aan;
dnS2O3dt = @(CS2O30,nS2O3,t,iAg,iAu,iPd) CS2O30*Fin-nS2O3*Fout/V-(2/F)*iAg*Acat-(2/F)*iAu*Acat-(4/(2*F))*iPd*Acat;

%modelling method vars
t(1) = 0; %time initial
h = 1/60; %step size, basically does a minute of time


for iter = 1:1:(tfinal/h)
    %calculated vars
    K(iter) = gammaAg*nAgS2O32(iter)+gammaAu*nAuS2O32(1)+gammaPd*2*nPdS2O34(1)+gammaH*nH(1);%update this
    R(iter) = cursivel/(Acat*K(1));
    %nernst potentials
    ErevAn(iter) = -1.23 + (Rgas*T/(4*F))*log(aO2*(gamH*nH(1)/V)^4);
    ErevFe(iter) = -0.77 + (Rgas*T/F)*log(gamFe3*(nFe2(iter)/V)/(gamFe2*(nFe3(iter)/V)));
    ErevAg(iter) = 0.060113 + (Rgas*T/(F))*log(((gamS2O3*(nS2O3(iter)/V))^2)/(gamAg*nAgS2O32(iter)/V));
    ErevAu(iter) = 0.153 + (Rgas*T/(F))*log(((gamS2O3*(nS2O3(iter)/V))^2)/(gamAg*nAuS2O32(iter)/V));
    ErevPd(iter) = 0.0862 + (Rgas*T/(2*F))*log(((gamS2O3*(nS2O3(iter)/V))^4)/((gamAg*nPdS2O34(iter))/V));
    ErevH(iter) = 0 + (Rgas*T/F)*log((aH2^0.5)/(gamH*nH(iter)));
    %solving for current and overpotentials
    CurrentFunc = @(eta) Aan*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T))) + iFe0*(exp(alphaFe*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T))-exp(-(1-alphaFe)*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T))))+Acat*( iAg0*(exp(alphaAg*F*eta(2)/(Rgas*T))-exp(-(1-alphaAg)*F*eta(2)/(Rgas*T))) + iAu0*(exp(alphaAu*F*(ErevAg(iter)+eta(2)-ErevAu(iter))/(Rgas*T))-exp(-(1-alphaAu)*F*(ErevAg(iter)+eta(2)-ErevAu(iter))/(Rgas*T))) + iPd0*(exp(alphaPd*2*F*(ErevAg(iter)+eta(2)-ErevPd(iter))/(Rgas*T))-exp(-(1-alphaPd)*2*F*(ErevAg(iter)+eta(2)-ErevPd(iter))/(Rgas*T)))+iH0*(exp(alphaH*F*(ErevAg(iter)+eta(2)-ErevH(iter))/(Rgas*T))-exp(-(1-alphaH)*F*(ErevAg(iter)+eta(2)-ErevH(iter))/(Rgas*T))));
    VappFunc = @(eta) abs(ErevAg(iter)-ErevAn(iter))+Aan*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T)))+ iFe0*(exp(alphaFe*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T))-exp(-(1-alphaFe)*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T))))*(R(1)+Rhardware)+abs(eta(1))+abs(eta(2))-Vapp;
    overpotentialFunc = @(eta) [CurrentFunc(eta);VappFunc(eta)];  
    etaGuess0 = [etaAn(iter) etaAg(iter)];
    options = optimoptions('fsolve','MaxFunctionEvaluations',1000,'MaxIterations',1000);
    eta = fsolve(overpotentialFunc,etaGuess0,options);
    %assigning overpotentials    
    etaAn(iter) = eta(1);
    etaFe(iter) = -ErevAn(iter)+eta(1)+ErevFe(iter);
    etaAg(iter) = eta(2);
    etaAu(iter) = ErevAg(iter)+eta(2)-ErevAu(iter);
    etaPd(iter) = ErevAg(iter)+eta(2)-ErevPd(iter);
    etaH(iter) = ErevAg(iter)+eta(2)-ErevH(iter);
    %assigning time and currents
    t(iter) = h*iter;
    iAg(iter) = iAg0*(exp(alphaAg*1*F*(etaAg(iter))/(Rgas*T))-exp(-(1-alphaAg)*1*F*etaAg(iter)/(Rgas*T)));
    iAu(iter) = iAu0*(exp(alphaAu*1*F*(etaAu(iter))/(Rgas*T))-exp(-(1-alphaAu)*1*F*etaAu(iter)/(Rgas*T)));
    iPd(iter) = iPd0*(exp(alphaPd*2*F*(etaPd(iter))/(Rgas*T))-exp(-(1-alphaPd)*2*F*etaPd(iter)/(Rgas*T)));
    iH(iter) = iH0*(exp(alphaH*1*F*(etaH(iter))/(Rgas*T))-exp(-(1-alphaH)*1*F*etaH(iter)/(Rgas*T)));
    iFe2(iter) = iFe0*(exp(alphaFe*F*(etaFe(iter))/(Rgas*T))-exp(-(1-alphaFe)*F*etaFe(iter)/(Rgas*T)));
    iAn(iter) = iAn0*(exp(alphaAn*4*F*(etaAn(iter))/(Rgas*T))-exp(-(1-alphaAn)*4*F*etaAn(iter)/(Rgas*T)));
    Itot(iter) = Acat*(iAg(iter)+iAu(iter)+iPd(iter)+iH(iter));
    Iantot(iter) = Aan*(iAn(iter)+iFe2(iter));
    Vappdiff(iter) = VappFunc(eta);
    Idiff(iter) = CurrentFunc(eta);
    %Ag
    k1 = dnAgdt(CAgS2O320(iter),nAgS2O32(iter),t(iter),iAg(iter));
    k2 = dnAgdt((CAgS2O320(iter)+(1/2)*h*k1),(nAgS2O32(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k1));
    k3 = dnAgdt((CAgS2O320(iter)+(1/2)*h*k2),(nAgS2O32(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k2));
    k4 = dnAgdt((CAgS2O320(iter)+h*k3),(nAgS2O32(iter)+h*k3),(t(iter)+h),(iAg(iter)+h*k3));
    nAgS2O32(iter+1) = subplus(nAgS2O32(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    CAgS2O320(iter+1) = CAgS2O320(iter);%nAg(iter+1);%for continuous recycle, use the commented out section. May add functionality for further coupled interactions.
    %Au
    k1 = dnAudt(CAuS2O320(iter),nAuS2O32(iter),t(iter),iAu(iter));
    k2 = dnAudt((CAuS2O320(iter)+(1/2)*h*k1),(nAuS2O32(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAu(iter)+(1/2)*h*k1));
    k3 = dnAudt((CAuS2O320(iter)+(1/2)*h*k2),(nAuS2O32(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAu(iter)+(1/2)*h*k2));
    k4 = dnAudt((CAuS2O320(iter)+h*k3),(nAuS2O32(iter)+h*k3),(t(iter)+h),(iAu(iter)+h*k3));
    nAuS2O32(iter+1) = subplus(nAuS2O32(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    CAuS2O320(iter+1) = CAuS2O320(iter);%nAu(iter+1);
    %Pd
    k1 = dnPddt(CPdS2O340(iter),nPdS2O34(iter),t(iter),iPd(iter));
    k2 = dnPddt((CPdS2O340(iter)+(1/2)*h*k1),(nPdS2O34(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iPd(iter)+(1/2)*h*k1));
    k3 = dnPddt((CPdS2O340(iter)+(1/2)*h*k2),(nPdS2O34(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iPd(iter)+(1/2)*h*k2));
    k4 = dnPddt((CPdS2O340(iter)+h*k3),(nPdS2O34(iter)+h*k3),(t(iter)+h),(iPd(iter)+h*k3));
    nPdS2O34(iter+1) = subplus(nPdS2O34(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    CPdS2O340(iter+1) = CPdS2O340(iter);%nPd(iter+1);
    %Fe2
    k1 = dnFe2dt(CFe20(iter),nFe2(iter),t(iter),iFe2(iter));
    k2 = dnFe2dt((CFe20(iter)+(1/2)*h*k1),(nFe2(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k1));
    k3 = dnFe2dt((CFe20(iter)+(1/2)*h*k2),(nFe2(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k2));
    k4 = dnFe2dt((CFe20(iter)+h*k3),(nFe2(iter)+h*k3),(t(iter)+h),(iFe2(iter)+h*k3));
    nFe2(iter+1) = subplus(nFe2(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    CFe20(iter+1) = CFe20(iter);%nPd(iter+1);
    %Fe3
    k1 = dnFe3dt(CFe30(iter),nFe3(iter),t(iter),iFe2(iter));
    k2 = dnFe3dt((CFe30(iter)+(1/2)*h*k1),(nFe3(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k1));
    k3 = dnFe3dt((CFe30(iter)+(1/2)*h*k2),(nFe3(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k2));
    k4 = dnFe3dt((CFe30(iter)+h*k3),(nFe3(iter)+h*k3),(t(iter)+h),(iFe2(iter)+h*k3));
    nFe3(iter+1) = subplus(nFe3(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    CFe30(iter+1) = CFe30(iter);%nPd(iter+1);
    %S2O3
    k1 = dnS2O3dt(CS2O30(iter),nS2O3(iter),t(iter),iAg(iter),iAu(iter),iPd(iter));
    k2 = dnS2O3dt((CS2O30(iter)+(1/2)*h*k1),(nS2O3(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k1),(iAu(iter)+(1/2)*h*k1),(iPd(iter)+(1/2)*h*k1));
    k3 = dnS2O3dt((CS2O30(iter)+(1/2)*h*k2),(nS2O3(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k2),(iAu(iter)+(1/2)*h*k2),(iPd(iter)+(1/2)*h*k2));
    k4 = dnS2O3dt((CS2O30(iter)+h*k3),(nS2O3(iter)+h*k3),(t(iter)+h),(iAg(iter)+h*k3),(iAu(iter)+h*k3),(iPd(iter)+h*k3));
    nS2O3(iter+1) = subplus(nS2O3(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    CS2O30(iter+1) = CS2O30(iter);%nPd(iter+1);
    %H
    nH(iter+1) = nH(iter); % refine this, im assuming something constantly balances pH here (or the conc is so big it dont matta)
    %other variables to calc
    etaAn(iter+1) = eta(1);
    etaAg(iter+1) = eta(2);
    %additional processing equations
    wAg(iter+1) = -(Acat*mmAg/F)*h*60*(sum(iAg)-iAg(end));
    wAu(iter+1) = -(Acat*mmAu/F)*h*60*(sum(iAu)-iAu(end));
    wPd(iter+1) = -(Acat*mmPd/(2*F))*h*60*(sum(iPd)-iPd(end));
end
t(iter+1)=h*(iter+1);
CAg = nAgS2O32./V;
CAu = nAuS2O32./V;
CPd = nPdS2O34./V;
RemAg = CAg./CAg(1);
RemAu = CAu./CAu(1);
RemPd = CPd./CPd(1);
subplot(2,2,1)
plot(t(1:end-1),Vappdiff)
title('Vapp Error')
subplot(2,2,2)
plot(t(1:end-1),Idiff)
title('Current error')
subplot(2,2,3)
plot(t,RemAg,t,RemAu,t,RemPd)
xlabel('Time (h)')
ylabel('% consumption from initial')
title('Remaining in solution over time, step size 1 min')
legend('Silver','Gold','Palladium')
subplot(2,2,4)
plot(t,wAg,t,wAu,t,wPd)
xlabel('Time (h)')
ylabel('Total amount deposited (g)')
title('Deposited solution over time, step size 1 min')
legend('Silver','Gold','Palladium')
toc

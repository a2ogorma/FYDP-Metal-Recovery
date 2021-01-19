clear all
tic
%%%initial cell concs%%%
%Corrosion
CAgS2O320(1) = 0.10;
CAuS2O320(1) = 0.05;
CPdS2O340(1) = 0.0001;
CH0(1) = 1e-10;
CFe20(1) = 0.005;
CFe30(1) = 0.001;
CS2O30(1) = 0.05;
COH0(1) = (1e-14)/CH0(1);
%Electrowinning
CAgS2O32(1) = 0.10;
CAuS2O32(1) = 0.05;
CPdS2O34(1) = 0.0001;
CH(1) = 1e-10;
CFe2(1) = 0.005;
CFe3(1) = 0.001;
CS2O3(1) = 0.05;
COH(1) = (1e-14)/CH0(1);

%%%Properties%%%
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

%%%System parameters%%%
Fin = 10; %L/s
Fout = Fin; %kept same for now, may want to make into different values to account for accumulation and controls
tfinal = 72;
%corrosion
Ecorr(1) = 0; %initial guess for corrosion potential
T_corr = 298.15;
V_corr = 100; %L
cursivel_corr = 10; %m, characteristic distance
A_corr = 10;%m2, area of uniform corrosion
Rhardware_corr = 0; %set to 0 for now before i do something with it
%electrowinning
Vapp = 8; %V
T_elec = 298.15; %K, we can play around with this but if we want to vary this kinetically then I oop
V_elec = 100; %L
cursivel_elec = 10; %m, characteristic distance
Acat_elec = 10; %m2, area
Aan_elec = 10; %m2, area
Rhardware_elec = 0; %set to 0 for now before i do something with it

%setting up initial conc in cell
nAgS2O32(1) = CAgS2O320(1)*V_elec; % moles of silver remaining
nAuS2O32(1) = CAuS2O320(1)*V_elec; % moles
nPdS2O34(1) = CPdS2O340(1)*V_elec; % moles
nS2O3(1) = CS2O30(1)*V_elec;
nFe2(1) = CFe20(1)*V_elec;
nFe3(1) = CFe30(1)*V_elec;
nH(1) = CH0(1)*V_elec; % moles
nOH(1) = COH0(1)*V_elec;
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
dnAgdt = @(CAg0,CAg,t,iAg) CAg0*Fin-CAg*Fout+(1/F)*iAg*Acat_elec;
dnAudt = @(CAu0,CAu,t,iAu) CAu0*Fin-CAu*Fout+(1/F)*iAu*Acat_elec;
dnPddt = @(CPd0,CPd,t,iPd) CPd0*Fin-CPd*Fout+(1/(2*F))*iPd*Acat_elec;
dnFe2dt = @(CFe20,CFe2,t,iFe2) CFe20*Fin-CFe2*Fout+(1/F)*iFe2*Aan_elec;
dnFe3dt = @(CFe30,CFe3,t,iFe2) CFe30*Fin-CFe3*Fout-(1/F)*iFe2*Aan_elec;
dnS2O3dt = @(CS2O30,CS2O3,t,iAg,iAu,iPd) CS2O30*Fin-CS2O3*Fout-(2/F)*iAg*Acat_elec-(2/F)*iAu*Acat_elec-(4/(2*F))*iPd*Acat_elec;

%modelling method vars
t(1) = 0; %time initial
h = 1/60; %step size, basically does a minute of time


for iter = 1:1:(tfinal/h)
    t(iter) = h*iter;
    %%%Corrosion%%%
    ErevAn_corr(iter) = 1.23 - (Rgas*T_corr/(4*F))*log(aO2*(gamH*CH0(iter))^4);
    ErevFe_corr(iter) = 0.77 - (Rgas*T_corr/F)*log(gamFe3*(CFe20(iter))/(gamFe2*(CFe30(iter))));
    ErevAg_corr(iter) = 0.060113 + (Rgas*T_corr/(F))*log(((gamS2O3*(CS2O30(iter)))^2)/(gamAg*CAgS2O320(iter)));
    ErevAu_corr(iter) = 0.153 + (Rgas*T_corr/(F))*log(((gamS2O3*(CS2O30(iter)))^2)/(gamAg*CAuS2O320(iter)));
    ErevPd_corr(iter) = 0.0862 + (Rgas*T_corr/(2*F))*log(((gamS2O3*(CS2O30(iter)))^4)/((gamAg*CPdS2O340(iter))));
    ErevH_corr(iter) = 0 + (Rgas*T_corr/F)*log((aH2^0.5)/(gamH*CH0(iter)));
    
    CorrosionFunc = @(Ecorr) A_corr*(i_BV((Ecorr-ErevFe_corr(iter)), iFe0, alphaFe, 1, T_corr) + i_BV((Ecorr-ErevAg_corr(iter)), iAg0, alphaAg, 1, T_corr) + i_BV((Ecorr-ErevAu_corr(iter)), iAu0, alphaAu, 1, T_corr) + i_BV((Ecorr-ErevPd_corr(iter)), iPd0, alphaPd, 2, T_corr) );
    Ecorr(iter) = fzero(CorrosionFunc,Ecorr(iter));
    Ecorr(iter+1) = Ecorr(iter); %setting up next initial guess
    iAg_corr(iter) = i_BV((Ecorr(iter)-ErevAg_corr(iter)), iAg0, alphaAg, 1, T_corr);
    iAu_corr(iter) = i_BV((Ecorr(iter)-ErevAu_corr(iter)), iAu0, alphaAu, 1, T_corr);
    iPd_corr(iter) = i_BV((Ecorr(iter)-ErevPd_corr(iter)), iPd0, alphaPd, 2, T_corr);
    iFe2_corr(iter) = i_BV((Ecorr(iter)-ErevFe_corr(iter)), iFe0, alphaFe, 1, T_corr);
    
    %Ag
    k1 = dnAgdt(CAgS2O32(iter),CAgS2O320(iter),t(iter),iAg_corr(iter));
    k2 = dnAgdt((CAgS2O32(iter)+(1/2)*h*k1),(CAgS2O320(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAg_corr(iter)+(1/2)*h*k1));
    k3 = dnAgdt((CAgS2O32(iter)+(1/2)*h*k2),(CAgS2O320(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAg_corr(iter)+(1/2)*h*k2));
    k4 = dnAgdt((CAgS2O32(iter)+h*k3),(CAgS2O320(iter)+h*k3),(t(iter)+h),(iAg_corr(iter)+h*k3));
    CAgS2O320(iter+1) = subplus(CAgS2O320(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Au
    k1 = dnAudt(CAuS2O32(iter),CAuS2O320(iter),t(iter),iAu_corr(iter));
    k2 = dnAudt((CAuS2O32(iter)+(1/2)*h*k1),(CAuS2O320(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAu_corr(iter)+(1/2)*h*k1));
    k3 = dnAudt((CAuS2O32(iter)+(1/2)*h*k2),(CAuS2O320(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAu_corr(iter)+(1/2)*h*k2));
    k4 = dnAudt((CAuS2O32(iter)+h*k3),(CAuS2O320(iter)+h*k3),(t(iter)+h),(iAu_corr(iter)+h*k3));
    CAuS2O320(iter+1) = subplus(CAuS2O320(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Pd
    k1 = dnPddt(CPdS2O34(iter),CPdS2O340(iter),t(iter),iPd_corr(iter));
    k2 = dnPddt((CPdS2O34(iter)+(1/2)*h*k1),(CPdS2O340(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iPd_corr(iter)+(1/2)*h*k1));
    k3 = dnPddt((CPdS2O34(iter)+(1/2)*h*k2),(CPdS2O340(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iPd_corr(iter)+(1/2)*h*k2));
    k4 = dnPddt((CPdS2O34(iter)+h*k3),(CPdS2O340(iter)+h*k3),(t(iter)+h),(iPd_corr(iter)+h*k3));
    CPdS2O340(iter+1) = subplus(CPdS2O340(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Fe2
    k1 = dnFe2dt(CFe2(iter),CFe20(iter),t(iter),iFe2_corr(iter));
    k2 = dnFe2dt((CFe2(iter)+(1/2)*h*k1),(CFe20(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iFe2_corr(iter)+(1/2)*h*k1));
    k3 = dnFe2dt((CFe2(iter)+(1/2)*h*k2),(CFe20(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iFe2_corr(iter)+(1/2)*h*k2));
    k4 = dnFe2dt((CFe2(iter)+h*k3),(CFe20(iter)+h*k3),(t(iter)+h),(iFe2_corr(iter)+h*k3));
    CFe20(iter+1) = subplus(CFe20(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Fe3
    k1 = dnFe3dt(CFe3(iter),CFe30(iter),t(iter),iFe2_corr(iter));
    k2 = dnFe3dt((CFe3(iter)+(1/2)*h*k1),(CFe30(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iFe2_corr(iter)+(1/2)*h*k1));
    k3 = dnFe3dt((CFe3(iter)+(1/2)*h*k2),(CFe30(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iFe2_corr(iter)+(1/2)*h*k2));
    k4 = dnFe3dt((CFe3(iter)+h*k3),(CFe30(iter)+h*k3),(t(iter)+h),(iFe2_corr(iter)+h*k3));
    CFe30(iter+1) = subplus(CFe30(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %S2O3
    k1 = dnS2O3dt(CS2O3(iter),CS2O30(iter),t(iter),iAg_corr(iter),iAu_corr(iter),iPd_corr(iter));
    k2 = dnS2O3dt((CS2O3(iter)+(1/2)*h*k1),(CS2O30(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAg_corr(iter)+(1/2)*h*k1),(iAu_corr(iter)+(1/2)*h*k1),(iPd_corr(iter)+(1/2)*h*k1));
    k3 = dnS2O3dt((CS2O3(iter)+(1/2)*h*k2),(CS2O30(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAg_corr(iter)+(1/2)*h*k2),(iAu_corr(iter)+(1/2)*h*k2),(iPd_corr(iter)+(1/2)*h*k2));
    k4 = dnS2O3dt((CS2O3(iter)+h*k3),(CS2O30(iter)+h*k3),(t(iter)+h),(iAg_corr(iter)+h*k3),(iAu_corr(iter)+h*k3),(iPd_corr(iter)+h*k3));
    CS2O30(iter+1) = subplus(CS2O30(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %H
    CH0(iter+1) = CH0(iter); % refine this, im assuming something constantly balances pH here (or the conc is so big it dont matta)
    
    %%%Electrowinning%%%
    %calculated vars
    K(iter) = gammaAg*CAgS2O32(iter)+gammaAu*CAuS2O32(iter)+gammaPd*2*CPdS2O34(iter)+gammaH*CH(iter);%update this
    R(iter) = cursivel_elec/(Acat_elec*K(1));
    %nernst potentials
    ErevAn(iter) = -1.23 + (Rgas*T_elec/(4*F))*log(aO2*(gamH*CH(iter))^4);
    ErevFe(iter) = -0.77 + (Rgas*T_elec/F)*log(gamFe3*(CFe2(iter))/(gamFe2*(CFe3(iter))));
    ErevAg(iter) = 0.060113 + (Rgas*T_elec/(F))*log(((gamS2O3*(CS2O3(iter)))^2)/(gamAg*CAgS2O32(iter)));
    ErevAu(iter) = 0.153 + (Rgas*T_elec/(F))*log(((gamS2O3*(CS2O3(iter)))^2)/(gamAg*CAuS2O32(iter)));
    ErevPd(iter) = 0.0862 + (Rgas*T_elec/(2*F))*log(((gamS2O3*(CS2O3(iter)))^4)/((gamAg*CPdS2O34(iter))));
    ErevH(iter) = 0 + (Rgas*T_elec/F)*log((aH2^0.5)/(gamH*CH(iter)));
    %solving for current and overpotentials
    CurrentFunc = @(eta) Aan_elec*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T_elec))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T_elec))) + iFe0*(exp(alphaFe*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T_elec))-exp(-(1-alphaFe)*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T_elec))))+Acat_elec*( iAg0*(exp(alphaAg*F*eta(2)/(Rgas*T_elec))-exp(-(1-alphaAg)*F*eta(2)/(Rgas*T_elec))) + iAu0*(exp(alphaAu*F*(ErevAg(iter)+eta(2)-ErevAu(iter))/(Rgas*T_elec))-exp(-(1-alphaAu)*F*(ErevAg(iter)+eta(2)-ErevAu(iter))/(Rgas*T_elec))) + iPd0*(exp(alphaPd*2*F*(ErevAg(iter)+eta(2)-ErevPd(iter))/(Rgas*T_elec))-exp(-(1-alphaPd)*2*F*(ErevAg(iter)+eta(2)-ErevPd(iter))/(Rgas*T_elec)))+iH0*(exp(alphaH*F*(ErevAg(iter)+eta(2)-ErevH(iter))/(Rgas*T_elec))-exp(-(1-alphaH)*F*(ErevAg(iter)+eta(2)-ErevH(iter))/(Rgas*T_elec))));
    VappFunc = @(eta) abs(ErevAg(iter)-ErevAn(iter))+Aan_elec*( iAn0*(exp(alphaAn*4*F*eta(1)/(Rgas*T_elec))-exp(-(1-alphaAn)*4*F*eta(1)/(Rgas*T_elec)))+ iFe0*(exp(alphaFe*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T_elec))-exp(-(1-alphaFe)*F*(-ErevAn(iter)+eta(1)+ErevFe(iter))/(Rgas*T_elec))))*(R(1)+Rhardware_elec)+abs(eta(1))+abs(eta(2))-Vapp;
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
    %assigning currents
    iAg(iter) = iAg0*(exp(alphaAg*1*F*(etaAg(iter))/(Rgas*T_elec))-exp(-(1-alphaAg)*1*F*etaAg(iter)/(Rgas*T_elec)));
    iAu(iter) = iAu0*(exp(alphaAu*1*F*(etaAu(iter))/(Rgas*T_elec))-exp(-(1-alphaAu)*1*F*etaAu(iter)/(Rgas*T_elec)));
    iPd(iter) = iPd0*(exp(alphaPd*2*F*(etaPd(iter))/(Rgas*T_elec))-exp(-(1-alphaPd)*2*F*etaPd(iter)/(Rgas*T_elec)));
    iH(iter) = iH0*(exp(alphaH*1*F*(etaH(iter))/(Rgas*T_elec))-exp(-(1-alphaH)*1*F*etaH(iter)/(Rgas*T_elec)));
    iFe2(iter) = iFe0*(exp(alphaFe*F*(etaFe(iter))/(Rgas*T_elec))-exp(-(1-alphaFe)*F*etaFe(iter)/(Rgas*T_elec)));
    iAn(iter) = iAn0*(exp(alphaAn*4*F*(etaAn(iter))/(Rgas*T_elec))-exp(-(1-alphaAn)*4*F*etaAn(iter)/(Rgas*T_elec)));
    Itot(iter) = Acat_elec*(iAg(iter)+iAu(iter)+iPd(iter)+iH(iter));
    Iantot(iter) = Aan_elec*(iAn(iter)+iFe2(iter));
    Vappdiff(iter) = VappFunc(eta);
    Idiff(iter) = CurrentFunc(eta);
    %Ag
    k1 = dnAgdt(CAgS2O320(iter),CAgS2O32(iter),t(iter),iAg(iter));
    k2 = dnAgdt((CAgS2O320(iter)+(1/2)*h*k1),(CAgS2O32(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k1));
    k3 = dnAgdt((CAgS2O320(iter)+(1/2)*h*k2),(CAgS2O32(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k2));
    k4 = dnAgdt((CAgS2O320(iter)+h*k3),(CAgS2O32(iter)+h*k3),(t(iter)+h),(iAg(iter)+h*k3));
    CAgS2O32(iter+1) = subplus(CAgS2O32(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Au
    k1 = dnAudt(CAuS2O320(iter),CAuS2O32(iter),t(iter),iAu(iter));
    k2 = dnAudt((CAuS2O320(iter)+(1/2)*h*k1),(CAuS2O32(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAu(iter)+(1/2)*h*k1));
    k3 = dnAudt((CAuS2O320(iter)+(1/2)*h*k2),(CAuS2O32(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAu(iter)+(1/2)*h*k2));
    k4 = dnAudt((CAuS2O320(iter)+h*k3),(CAuS2O32(iter)+h*k3),(t(iter)+h),(iAu(iter)+h*k3));
    CAuS2O32(iter+1) = subplus(CAuS2O32(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Pd
    k1 = dnPddt(CPdS2O340(iter),CPdS2O34(iter),t(iter),iPd(iter));
    k2 = dnPddt((CPdS2O340(iter)+(1/2)*h*k1),(CPdS2O34(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iPd(iter)+(1/2)*h*k1));
    k3 = dnPddt((CPdS2O340(iter)+(1/2)*h*k2),(CPdS2O34(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iPd(iter)+(1/2)*h*k2));
    k4 = dnPddt((CPdS2O340(iter)+h*k3),(CPdS2O34(iter)+h*k3),(t(iter)+h),(iPd(iter)+h*k3));
    CPdS2O34(iter+1) = subplus(CPdS2O34(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Fe2
    k1 = dnFe2dt(CFe20(iter),CFe2(iter),t(iter),iFe2(iter));
    k2 = dnFe2dt((CFe20(iter)+(1/2)*h*k1),(CFe2(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k1));
    k3 = dnFe2dt((CFe20(iter)+(1/2)*h*k2),(CFe2(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k2));
    k4 = dnFe2dt((CFe20(iter)+h*k3),(CFe2(iter)+h*k3),(t(iter)+h),(iFe2(iter)+h*k3));
    CFe2(iter+1) = subplus(CFe2(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %Fe3
    k1 = dnFe3dt(CFe30(iter),CFe3(iter),t(iter),iFe2(iter));
    k2 = dnFe3dt((CFe30(iter)+(1/2)*h*k1),(CFe3(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k1));
    k3 = dnFe3dt((CFe30(iter)+(1/2)*h*k2),(CFe3(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iFe2(iter)+(1/2)*h*k2));
    k4 = dnFe3dt((CFe30(iter)+h*k3),(CFe3(iter)+h*k3),(t(iter)+h),(iFe2(iter)+h*k3));
    CFe3(iter+1) = subplus(CFe3(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %S2O3
    k1 = dnS2O3dt(CS2O30(iter),CS2O3(iter),t(iter),iAg(iter),iAu(iter),iPd(iter));
    k2 = dnS2O3dt((CS2O30(iter)+(1/2)*h*k1),(CS2O3(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k1),(iAu(iter)+(1/2)*h*k1),(iPd(iter)+(1/2)*h*k1));
    k3 = dnS2O3dt((CS2O30(iter)+(1/2)*h*k2),(CS2O3(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k2),(iAu(iter)+(1/2)*h*k2),(iPd(iter)+(1/2)*h*k2));
    k4 = dnS2O3dt((CS2O30(iter)+h*k3),(CS2O3(iter)+h*k3),(t(iter)+h),(iAg(iter)+h*k3),(iAu(iter)+h*k3),(iPd(iter)+h*k3));
    CS2O3(iter+1) = subplus(CS2O3(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    %H
    CH(iter+1) = CH(iter); % refine this, im assuming something constantly balances pH here (or the conc is so big it dont matta)
    %other variables to calc
    etaAn(iter+1) = eta(1);
    etaAg(iter+1) = eta(2);
    %additional processing equations
    wAg(iter+1) = -(Acat_elec*mmAg/F)*h*60*(sum(iAg)-iAg(end));
    wAu(iter+1) = -(Acat_elec*mmAu/F)*h*60*(sum(iAu)-iAu(end));
    wPd(iter+1) = -(Acat_elec*mmPd/(2*F))*h*60*(sum(iPd)-iPd(end));
end
t(iter+1)=h*(iter+1);
CAg = nAgS2O32./V_elec;
CAu = nAuS2O32./V_elec;
CPd = nPdS2O34./V_elec;
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

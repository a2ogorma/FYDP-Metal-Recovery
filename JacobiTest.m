clear all
%inlet concs for now
CAg0(1) = 0.10;
CAu0(1) = 0.05;
CPd0(1) = 0.0001;
CH0(1) = 1e-10;

%conductivity values, check on this (S m2 /mol)
gammaAg = 61.90e-4;
gammaAu = 4.1e7;
gammaPd = 0.1;
gammaH = 349.81e-4;
%divergence from ideality values (check on this shit too but low
%priority)... something about fugacity comes to mind???
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
Vapp = 3; %V
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

%calculated vars (will be iterated)
K(1) = gammaAg*nAg(1)+gammaAu*nAu(1)+gammaPd*2*nPd(1)+gammaH*nH(1);
R(1) = cursivel/(Acat*K(1));
Ean(1) = -1.23 - (Rgas*T/(4*F))*log(aO2*(gamH*nH(1)/V)^4);
ErevAg(1) = 0.7989 + (Rgas*T/(F))*log(1/(gamAg*nAg(1)/V));
ErevAu(1) = 1.692 + (Rgas*T/(F))*log(1/(gamAg*nAu(1)/V));
ErevPd(1) = 0.951 + (Rgas*T/(2*F))*log(1/(((gamAg*nPd(1))/V)^2));


%assume well mixed, or bulk and surface concentrations are the same
alphaAn = 0.5;
alphaAg = 0.5; %base assumption as most systems are close to this number
alphaAu = 0.5; %base assumption as most systems are close to this number
alphaPd = 0.5; %base assumption as most systems are close to this number
iAg0 = 7.3e-9; %A/m2, from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1007/s10800-007-9434-x - verify this link applies
iAu0 = 2E-8; %A/m2, value of dissolution including Na2S in dissolution reaction from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1134/S1023193506040021, converting A/cm2 to A/m2 
iPd0 = 4e-9; %using a mix of a few from that reference
iAn0 = 1E-8; %A/cm2

%initial guesses
I(1) = 9000000;%abs((Vapp - abs(ErevAg(1)-Ean(1)) - 3)/(R(1)+Rhardware)); %this is the initial initial guess for I
etaAn(1) = 1;
etaAg(1) = -1;
etaAu(1) = -1;
etaPd(1) = -1;
etaAn(1) = Rgas*T*log((I(1)/Aan)/iAn0)/(4*F); %solving for anodic overpotential from estimated current
anodefunc = @(etaAn)Aan*iAn0*(exp(alphaAn*4*F*etaAn/(Rgas*T))-exp(-(1-alphaAn)*4*F*etaAn/(Rgas*T)))-I(1);
       
etaAn(1) = fzero(anodefunc,etaAn(1));
for iter = 1:1:500
    %function
    f(1) = abs(ErevAg(1)-Ean(1))+I(1)*(R(1)+Rhardware)+abs(etaAn(iter))+abs(etaAg(iter))-Vapp;
    f(2) = Aan*iAn0*(exp(alphaAn*4*F*etaAn(iter)/(Rgas*T))-exp(-(1-alphaAn)*4*F*etaAn(iter)/(Rgas*T)))-I(1);
    f(3) = Acat*(iAg0*(exp(alphaAg*F*etaAg(iter)/(Rgas*T))-exp(-(1-alphaAg)*F*etaAg(iter)/(Rgas*T)))+iAu0*(exp(alphaAu*F*etaAu(iter)/(Rgas*T))-exp(-(1-alphaAu)*F*etaAu(iter)/(Rgas*T)))+iPd0*(exp(alphaPd*2*F*etaPd(iter)/(Rgas*T))-exp(-(1-alphaPd)*2*F*etaPd(iter)/(Rgas*T))))+I(1);
    f(4) = ErevAu(1)+abs(etaAu(iter))-ErevAg(1)-abs(etaAg(iter));
    f(5) = ErevPd(1)+abs(etaPd(iter))-ErevAg(1)-abs(etaAg(iter));

    J(1,1) = (R(1)+Rhardware);
    J(1,2) = abs(etaAn(iter))/etaAn(iter);
    J(1,3) = abs(etaAg(iter))/etaAg(iter);
    J(1,4) = 0;
    J(1,5) = 0;
    J(2,1) = -1;
    J(2,2) = Aan*iAn0*((alphaAn*4*F/(Rgas*T))*exp(alphaAn*4*F*etaAn(iter)/(Rgas*T))+((1-alphaAn)*4*F/(Rgas*T))*exp(-(1-alphaAn)*4*F*etaAn(iter)/(Rgas*T)));
    J(2,3) = 0;
    J(2,4) = 0;
    J(2,5) = 0;
    J(3,1) = 1;
    J(3,2) = 0;
    J(3,3) = Acat*iAg0*((alphaAg*F/(Rgas*T))*exp(alphaAg*F*etaAg(iter)/(Rgas*T))+((1-alphaAg)*F/(Rgas*T))*exp(-(1-alphaAg)*F*etaAg(iter)/(Rgas*T)));
    J(3,4) = Acat*iAu0*((alphaAu*F/(Rgas*T))*exp(alphaAu*F*etaAu(iter)/(Rgas*T))+((1-alphaAu)*F/(Rgas*T))*exp(-(1-alphaAu)*F*etaAu(iter)/(Rgas*T)));
    J(3,5) = Acat*iPd0*((alphaPd*2*F/(Rgas*T))*exp(alphaPd*2*F*etaPd(iter)/(Rgas*T))+((1-alphaPd)*2*F/(Rgas*T))*exp(-(1-alphaPd)*2*F*etaPd(iter)/(Rgas*T)));
    J(4,1) = 0;
    J(4,2) = 0;
    J(4,3) = -1;
    J(4,4) = 1;
    J(4,5) = 0;
    J(5,1) = 0;
    J(5,2) = 0;
    J(5,3) = -1;
    J(5,4) = 0;
    J(5,5) = 1;

    delX = (-f')\J';
    I(iter+1) = I(iter)+delX(1);
    etaAn(iter+1) = etaAn(iter)+delX(2);
    etaAg(iter+1) = etaAg(iter)+delX(3);
    etaAu(iter+1) = etaAu(iter)+delX(4);
    etaPd(iter+1) = etaPd(iter)+delX(5);
end
disp(I(end))%abs((Vapp - abs(ErevAg(1)-Ean(1)) - 3)/(R(1)+Rhardware)); %this is the initial initial guess for I
disp(etaAn(end))
disp(etaAg(end))
disp(etaAu(end))
disp(etaPd(end))
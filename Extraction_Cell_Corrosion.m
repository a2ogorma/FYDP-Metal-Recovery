%{
	Continuous stirred tank model for base metal recovery cell. Corrosion
	Edition
	Author: Matt/Ryan
	To Add: Applied Voltage calc. Non-constant volume? Mass transfer stuff?
	Multiple rxns
%} 
clear all 
%universal constants
R = 8.314; %J/(mol K)
F = 96485.3329; %C/mol
% general parameters
temp = 273; %K
pressure = 1; % atm
volume = 10; %L
Q = 0.01; % L/s (flowrate)
%{
	Reactions
	Anodic: Fe2+ --> Fe3+ + e-
	Cathodic: Cu2+ + 2e- --> Cu(s)
%}
% # of electrons in rxn
z_Cu = 2;
z_Fe = 1;
% # of total electrons in rxn
n = 3;
% exchange current densities
i0_Cu = 2E-5; %A/cm^2
i0_Fe = 1E-8; %A/cm^2
% charge transfer coefficients
alpha_Cu = 0.5;
alpha_Fe = 0.5;
% Standard half reaction potentials vs. SHE @ 298 K, 1 atm, 1 M https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
Eo_Cu = 0.337; %V
Eo_Fe = 0.77; %V
%activity coefficients of ions in solution
gamma_Cu2 = 1;
gamma_Fe2 = 1;
gamma_Fe3 = 1;
%initial concentrations in mol/L
%Outlet (In CSTR)
Ci_Cu2 = 0.01;
Ci_Fe2 = 0.01; %Changed these to be able to make the Fe equations work
Ci_Fe3 = 0.01; 
%inlet concentrations
C_Cu2_o = 0.7;
C_Fe2_o = 0.4;
C_Fe3_o = 0.1;

%Electrode Potentials
E_cat = 0; %V
E_an = 1.5; %V

%Electrode areas
S_cat = 100; %cm^2
S_an = 100; %cm^2

%Specific Surface Area
V_bed = 1; %m3 Volume of bed holding the particles assuming the bed is completly full.
radius = 0.001; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
SSA = 3/radius; %m2/m3 Specific Surface area of spheres.
density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
SA = V_bed*density*SSA; %m2 Total surface area of spheres in bed.

%Mass balance on generic ion in reactor
dC_dt = @(C, C_o, E, Erev, a, z, S, i0)((C_o-C)*Q + i0*(exp((a*z*F)/(R*temp)*(E-Erev))-exp(-(a*z*F)/(R*temp)*(E-Erev))*S/F/z)/volume);

%Solver
tfinal = 600; %total time in seconds
h = 1; %Step size = 1 s

%initializing concentrations
C_Cu2 = zeros(1, tfinal/h+1);
C_Fe2 = zeros(1, tfinal/h+1);
C_Fe3 = zeros(1, tfinal/h+1);
x0 = zeros(1, tfinal/h+1);
C_Cu2(1) = Ci_Cu2;
C_Fe2(1) = Ci_Fe2;
C_Fe3(1) = Ci_Fe3;
i_Cu = zeros(1, tfinal/h+1);
i_Fe = zeros(1, tfinal/h+1);
i_Fe(1) = 0;
i_Cu(1) = 0; %Literally just put this in as placeholder. Spend a couple minutes looking over appropriate value.

%Nernst potential arrays
Erev_Cu = zeros(1, tfinal/h+1);
Erev_Fe = zeros(1, tfinal/h+1);
Ecorr = zeros(1, tfinal/h+1);
icorr = zeros(1, tfinal/h+1);
Erev_Cu(1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2(1)));
Erev_Fe(1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2(1)/gamma_Fe3*C_Fe3(1));

x = ((n/z_Fe)*i0_Fe*(SA)/((n/z_Cu)*i0_Cu*(SA))); %Combining this here to save space later
y_an = (alpha_Fe*z_Fe*F)/(R*temp);
y_cat = ((1-alpha_Cu)*z_Cu*F)/(R*temp);

Ecorr(1) = ((log(x)*y_an*Erev_Fe(1))+(y_cat*Erev_Cu(1)))/((log(x)*y_an)+y_cat);
icorr(1) = (n/z_Cu)*i0_Cu*(exp((1-alpha_Cu)*z_Cu*F*(Ecorr(1)-Erev_Cu(1))/(R*temp)));

for i = 1:1:(tfinal/h) 
	%Runge Kutta 4 for Cu
	k1 = dC_dt(C_Cu2(i), C_Cu2_o, E_cat, Erev_Cu(i), alpha_Cu, z_Cu, S_cat,i0_Cu);
	k2 = dC_dt(C_Cu2(i)+h*k1/2, C_Cu2_o, E_cat, Erev_Cu(i), alpha_Cu, z_Cu, S_cat,i0_Cu);
	k3 = dC_dt(C_Cu2(i)+h*k2/2, C_Cu2_o, E_cat, Erev_Cu(i), alpha_Cu, z_Cu, S_cat,i0_Cu);
	k4 = dC_dt(C_Cu2(i)+h*k3, C_Cu2_o, E_cat, Erev_Cu(i), alpha_Cu, z_Cu, S_cat,i0_Cu);
	C_Cu2(i+1) = C_Cu2(i) + 1/6*h*(k1+2*k2+2*k3+k4); %Runge Kutta for Mass balance
    %Runge Kutta 4 for Fe2
	k1 = dC_dt(C_Fe2(i), C_Fe2_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	k2 = dC_dt(C_Fe2(i)+h*k1/2, C_Fe2_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	k3 = dC_dt(C_Fe2(i)+h*k2/2, C_Fe2_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	k4 = dC_dt(C_Fe2(i)+h*k3, C_Fe2_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	C_Fe2(i+1) = C_Fe2(i) + 1/6*h*(k1+2*k2+2*k3+k4); %Runge Kutta for Mass balance
    %Runge Kutta 4 for Fe3
	k1 = dC_dt(C_Fe3(i), C_Fe3_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	k2 = dC_dt(C_Fe3(i)+h*k1/2, C_Fe3_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	k3 = dC_dt(C_Fe3(i)+h*k2/2, C_Fe3_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	k4 = dC_dt(C_Fe3(i)+h*k3, C_Fe3_o, E_an, Erev_Fe(i), alpha_Fe, z_Fe, S_an,i0_Fe);
	C_Fe3(i+1) = C_Fe3(i) + 1/6*h*(k1+2*k2+2*k3+k4); %Runge Kutta for Mass balance
    %i0_Cu*(exp((alpha_Cu*z_Cu*F)/(R*temp)*(E_cat-Erev_Cu(i)))-exp(-(alpha_Cu*z_Cu*F)/(R*temp)*(E_cat-Erev_Cu(i)))*S_cat/F/z_Cu)/volume; %Alternate BV Eq'n
    %Nernst
	Erev_Cu(i+1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2(i+1)));
	Erev_Fe(i+1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2(i+1)/gamma_Fe3*C_Fe3(i+1));
    
    %Ecorr
    fun = @(Ecorr)cor(Ecorr, [Erev_Cu(i+1), Erev_Fe(i+1)], [i0_Cu, i0_Fe], [alpha_Cu, alpha_Fe],[z_Cu, z_Fe], temp);
    Ecorr(i+1) = fzero(fun,x0(i));
    x0(i+1) = Ecorr(i+1);
    icorr(i+1) = i0_Cu*(z_Cu*F/(R*temp))*(Ecorr(i+1)-Erev_Cu(i+1));
    %{ 
    Old overpotential stuff that did not use BV
    %Ecorr(i+1) = ((log(x)*y_an*Erev_Fe(i+1))+(y_cat*Erev_Cu(i+1)))/((log(x)*y_an)+y_cat); %Overpotential |n| > 0.2.
    %icorr(i+1) = (n/z_Cu)*i0_Cu*(exp((1-alpha_Cu)*z_Cu*F*(Ecorr(i+1)-Erev_Cu(i+1))/(R*temp))); %Overpotential |n| > 0.2
    %Ecorr(i+1) = ((Erev_Fe(i+1)*(((i0_Fe*z_Fe)/(i0_Cu*z_Cu))-Erev_Cu(i+1)))/(((i0_Fe*z_Fe)/(i0_Cu*z_Cu))-1)); %Overpotential |n| < 0.2
    %icorr(i+1) = i0_Cu*(z_Cu*F/(R*temp))*(Ecorr(i+1)-Erev_Cu(i+1)); %Overpotential |n| < 0.2
    %}
end
i = (1:1:(tfinal/h+1)); %Didn't want to deal with i in the 'for' function. Defined it the same here.
figure
subplot(2,2,1)
    plot(i,Erev_Cu)
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    xlim([0 tfinal/h+1])
    title('Erev Cu over time')
    legend('Copper')
subplot(2,2,2)
    plot(i,Erev_Fe)
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    xlim([0 tfinal/h+1])
    title('Erev Fe over time')
    legend('Iron')
subplot(2,2,3)
    plot(i,Ecorr)
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    xlim([0 tfinal/h+1])
    title('Corrosion Voltage over time')
subplot(2,2,4)
    plot(i,icorr)
    xlabel('Time (s)')
    ylabel('Current (A)')
    xlim([0 tfinal/h+1])
    title('Corrosion Current over time')
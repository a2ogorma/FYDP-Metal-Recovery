%{
	Continuous stirred tank model for base metal extraction/recovery cell
    Using ode45
	To Add: Applied Voltage calc. Non-constant volume?
    Order of compounds: Inert Cu Sn Al Pb Fe Zn Ca Ni Ag Au Pd
%} 
clear variables

propertiesBaseMetals %calls all the property values as set in the applicable file. Note this is currently only applicable to base metals

% system parameters
temp = 298; %K
pres = 1; % atm
vol_cell = 250; %L
Q = 5;%; % L/s (flowrate)
%Electrode areas
S_cat = 500; %cm^2
S_an = 500; %cm^2
%Cross sectional area of cell
A_cell = 500; %cm^2
%Length b/w electrodes
l = 100; %cm
%Applied Voltage (potentiostat)
V_app = 12; %V
%Extraction vessel parameters
vol_lch = 200; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
m_PCB_total = 80; %kg Mass of crushed PCBs
r_particles = 0.001; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
tfinal = 0.10; %s

%Weight fraction composition of PCB
%Inert Cu Sn Al Pb Fe 
wtfrac_PCB(1:6) = [0.783 0.151 1.80E-2 1.57E-2 1.16E-2 0.781E-2]; 
%Ag Au Pd
wtfrac_PCB(7:9) = [0.0130E-2 0.00580E-2 0.00286E-2];
%masses of each component:
m_PCB = m_PCB_total*wtfrac_PCB;
%Convert to volume and volume fraction
V_PCB = m_PCB./rho; %Element volumes in m^3
vfrac_PCB = V_PCB/sum(V_PCB);
V_PCB_total = sum(V_PCB);

packing_density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
fill_pct = 100*sum(V_PCB)/(0.001*packing_density*vol_lch);
if fill_pct > 75
    txt = ['Warning: PCB mass is ', num2str(fill_pct), '% of total lching vessel volume'];
    disp(txt);
end
%Surface area calculation for corrosion
SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
S_PCB_total = V_PCB_total*packing_density*SSA; %m2 Total surface area of spheres in lch.
S_PCB = vfrac_PCB*S_PCB_total; %array with exposed surface area of each component
S_PCB = S_PCB*100^2; %convert m^2 to cm^2

%Calculate number of crushed PCB particles
n_particles = sum(V_PCB)*3/(4*pi*r_particles^3);
%initial concentrations in mol/L
%Cell Concentrations (recovery)
Ci_Cu2_cell = 0;
Ci_Sn2_cell = 0.0;
Ci_Al3_cell = 0.0;
Ci_Pb2_cell = 0.0;
Ci_Fe2_cell = 0.5;
Ci_Fe3_cell = 0.1;
Ci_Ag_cell = 0.0;
Ci_Au_cell = 0.0;
Ci_Pd2_cell = 0.0;
Ci_H_cell = 0.5;
Ci_Cl_cell = 2*(Ci_Cu2_cell+Ci_Fe2_cell)+Ci_H_cell+3*Ci_Fe3_cell; %Calculation to ensure 
%electrolyte has net neutral charge
Ci_cell = [Ci_Cu2_cell Ci_Sn2_cell Ci_Al3_cell Ci_Pb2_cell Ci_Fe2_cell ...
    Ci_Fe3_cell Ci_Ag_cell Ci_Au_cell Ci_Pd2_cell Ci_H_cell Ci_Cl_cell];

%leching vessel concentrations (extraction)
Ci_Cu2_lch = 0;
Ci_Sn2_lch = 0;
Ci_Al3_lch = 0;
Ci_Pb2_lch = 0;
Ci_Fe2_lch = 0.5;
Ci_Fe3_lch = 0.1;
Ci_Ag_lch = 0;
Ci_Au_lch = 0;
Ci_Pd2_lch = 0;
Ci_H_lch = 0.5;
Ci_Cl_lch = 2*(Ci_Cu2_lch+Ci_Fe2_lch)+Ci_H_lch+3*Ci_Fe3_lch;
Ci_lch = [Ci_Cu2_lch Ci_Sn2_lch Ci_Al3_lch Ci_Pb2_lch Ci_Fe2_lch ... 
    Ci_Fe3_lch Ci_Ag_lch Ci_Au_lch Ci_Pd2_lch Ci_H_lch Ci_Cl_lch];

%initializing solution concentration and solid mass vector
Cm_i = [Ci_cell Ci_lch m_PCB];

%solve concentration profiles
tspan = [0 tfinal];
options = odeset('NonNegative',1:31);
balance_solver = @(t, Cm) ion_balance(t, Cm, temp, pres, vol_cell, vol_lch, Q, S_an, S_cat, V_app, n_particles, l, A_cell);
[t, Cm] = ode15s(balance_solver, tspan, Cm_i);

for j = 1:1:length(t)
    disp(t(j))
    %Nernst Potentials
    CmStep = Cm(j,:);
    [Erev_cell(j,:), Erev_lch(j,:)] = nernstPotential(CmStep,temp,aH2,aO2);
    %resistance calculation for IR drops
    kappa = 1000*(CmStep(1)*lamda(1)+CmStep(2)*lamda(2)+CmStep(3)*lamda(3)+...
        CmStep(4)*lamda(4)+CmStep(5)*lamda(5)+CmStep(6)*lamda(6)+CmStep(7)*lamda(7)+...
        CmStep(8)*lamda(8)+CmStep(9)*lamda(9)+CmStep(10)*lamda(10)+CmStep(11)*lamda(11));
    r_sol(j) = l/A_cell/kappa*100; %ohms
    r_hardware = 1; %ohms
    
    %solve cell currents and electrode potentials
    onCathode = [1 1 1 1 0 1 1 1 1 1 0];
    onAnode = -(onCathode-1);
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol(j), r_hardware, Erev_cell(j,:), onCathode, S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [0.2, 0.1, -0.1];
    options = optimoptions(@fsolve, 'Display','final', 'MaxFunctionEvaluations', 3000);
    x = fsolve(solver, x0, options);
    I(j) = x(1);
    E_an(j) = x(2);
    E_cat(j) = x(3);
    eta_cat(j,:) = E_cat(j) - Erev_cell(j,:);
    eta_an(j,:) = E_an(j) - Erev_cell(j,:);
    i_cat(j,:) = onCathode.*i_BV(eta_cat(j,:), i0, alphas, z, temp);
    I_cat(j,:) = i_cat(j,:)*S_cat;
    i_an(j,:) = onAnode.*i_BV(eta_an(j,:), i0, alphas, z, temp);
    I_an(j,:) = i_an(j,:)*S_an;
    I_cell(j,:) = I_cat(j,:)+I_an(j,:); %overall current for rxn i in cell
    
    m_PCB(j,:) = (Cm(j,23:31));
    m_PCB_total(j) = sum(m_PCB(j,:));
    wtfrac_PCB(j,:) = m_PCB(j,:)/m_PCB_total(j);
    %Convert weight to volume and volume fraction
    V_PCB(j,:) = m_PCB(j,:)./rho; %Element volumes in m^3
    vfrac_PCB(j,:) = V_PCB(j,:)/sum(V_PCB(j,:));
    V_PCB_total(j) = sum(V_PCB(j,:));
    
    %calculate particle radius in m
    r_particles(j) = (V_PCB_total(j)/(n_particles*4*pi)*3)^(1/3);
    
    %Surface area calculation for corrosion
    packing_density = 0.6;
    SSA(j) = 3/r_particles(j); %m2/m3 Specific Surface area of spheres.
    S_PCB_total(j) = V_PCB_total(j)*packing_density*SSA(j); %m2 Total surface area of spheres in lch.
    S_PCB(j,:) = vfrac_PCB(j,:)*S_PCB_total(j); %array with exposed surface area of each component
    S_PCB(j,:) = subplus(S_PCB(j,:)*100^2); %convert m^2 to cm^2
    
    %solve extraction bed corrosion rate
    j0 = 0; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, Erev_lch(j,:), S_PCB(j,:), temp);
    E_corr(j) = fzero(cor_solver, j0);
    i_corr(j,:) = i_BV(E_corr(j)-Erev_lch(j,:), i0, alphas, z, temp);
    S_corr(j,:) = [S_PCB(j,2:5) sum(S_PCB(j,:)) S_PCB(j,6:9) 0 sum(S_PCB(j,:))];
    I_corr(j,:) = S_corr(j,:).*i_corr(j,:);
end
%{
%Calculate currents, overpotentials
eta_Cu = E_cat - Erev_Cu_cell;
eta_Sn = E_cat - Erev_Sn_cell;
eta_Ni = E_cat - Erev_Ni_cell;
I_Cu = S_cat*i_BV(eta_Cu, i0_Cu, alpha_Cu, z_Cu, temp);
I_Sn = S_cat*i_BV(eta_Sn, i0_Sn, alpha_Sn, z_Sn, temp);
I_Ni = S_cat*i_BV(eta_Cu, i0_Ni, alpha_Ni, z_Ni, temp);

%plots
subplot(3,2,1)
plot(t,Cm(:,1),t,Cm(:,8))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Copper')
subplot(3,2,2)
plot(t,Cm(:,2),t,Cm(:,9))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(II)')
subplot(3,2,3)
plot(t,Cm(:,3),t,Cm(:,10))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Iron(III)')
subplot(3,2,4)
plot(t,Cm(:,4),t,Cm(:,11))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Tin')
subplot(3,2,5)
plot(t,Cm(:,5),t,Cm(:,12))
legend('Cell','Bed')
xlabel('Time (s)')
ylabel('Concentrations (M)')
title('Nickel')
subplot(3,2,6)
plot(t,I_corr,t,I_Cu,t,I_Fe)
xlabel('Time (s)')
ylabel('Current (A)')
title('Currents')
legend('Corr','CopperCell','IronCell')
%}

function func = cor(Ecorr, Erev, S_PCB, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alphas z
    i_corr = i_BV(Ecorr-Erev, i0, alphas, z, temp);
    %arrange surface areas in proper order
    S_corr = [S_PCB(2:5) sum(S_PCB) S_PCB(6:9) 0 sum(S_PCB)];
    if abs(sum(S_corr.*i_corr)) == Inf
        disp(S_corr)
        disp(i_corr)
        error("Corrosion currents infinite");
    end
    %imag(sum(S_corr.*i_corr))
    func = sum(S_corr.*i_corr);
end
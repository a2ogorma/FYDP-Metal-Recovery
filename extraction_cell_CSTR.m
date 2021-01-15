%{
	Continuous stirred tank model for base metal recovery cell
	Author: Matt/Ryan
	To Add: Applied Voltage calc. Non-constant volume?
%} 
clear all
global R F;
%universal constants
R = 8.314; %J/(mol K)
F = 96485.3329; %C/mol
% system parameters
global temp;
temp = 298; %K
pressure = 1; % atm
volume = 400; %L
Q = 1; % L/s (flowrate)
%{
	Reactions
	Anodic: Fe2+ --> Fe3+ + e-
	Cathodic: Cu2+ + 2e- --> Cu(s)
%}
%Defining constants as global variables for use in cell solver function
global z_Cu z_Fe i0_Cu i0_Fe alpha_Cu alpha_Fe
% # of electrons in rxn
z_Cu = 2;
z_Fe = 1;
% exchange current densities
i0_Cu = 5E-5; %A/cm^2
i0_Fe = 5E-8; %A/cm^2
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
Ci_Cu2 = 0.1;
Ci_Fe2 = 0.1;
Ci_Fe3 = 0.001; 
%inlet concentrations
C_Cu2_o = 0.7;
C_Fe2_o = 0.7;
C_Fe3_o = 0.05;

%Electrode areas
S_cat = 100; %cm^2
S_an = 100; %cm^2
%Applied Voltage (potentiostat)
V_app = 4; %V

%Mass balance on generic ion in reactor
dC_dt = @(C, C_o, I, z) ((C_o-C)*Q + I*F/z)/volume;

%Solver
tfinal = 0.001; %total time in seconds
h = 0.00001; %Step size = 1 s
t = [0:h:tfinal];
%initializing concentrations
C_Cu2 = zeros(1, tfinal/h+1);
C_Fe2 = zeros(1, tfinal/h+1);
C_Fe3 = zeros(1, tfinal/h+1);
C_Cu2(1) = Ci_Cu2;
C_Fe2(1) = Ci_Fe2;
C_Fe3(1) = Ci_Fe3;
%Nernst potential arrays
Erev_Cu = zeros(1, tfinal/h+1);
Erev_Fe = zeros(1, tfinal/h+1);
Erev_Cu(1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2(1)));
Erev_Fe(1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2(1)/gamma_Fe3*C_Fe3(1));
%Resistances
r_sol = ones(1, tfinal/h); %ohms, add calculation for this later
r_hardware = 0; %ohms
%Currents, potential, overpotential arrays
I_Cu = zeros(1, tfinal/h+1);
I_Fe = zeros(1, tfinal/h+1);
E_an = zeros(1, tfinal/h+1);
E_cat = zeros(1, tfinal/h+1);
for k = 1:1:(tfinal/h)
    %solve for currents and electrode potentials using fsolve
    %wrapper function for cell solver
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol(k), r_hardware, Erev_Fe(k), Erev_Cu(k), S_an, S_cat);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    x = fsolve(solver, x0);
    I_an = x(1);
    E_an(k) = x(2);
    E_cat(k) = x(3);
    I_Cu(k) = - I_an;
    I_Fe(k) = I_an;
    %4th order runge-kutta to solve mass balance differential equations
	k1 = dC_dt(C_Cu2(k), C_Cu2_o, I_Cu(k), z_Cu);
	k2 = dC_dt(C_Cu2(k)+h*k1/2, C_Cu2_o, I_Cu(k), z_Cu);
	k3 = dC_dt(C_Cu2(k)+h*k2/2, C_Cu2_o, I_Cu(k), z_Cu);
	k4 = dC_dt(C_Cu2(k)+h*k3, C_Cu2_o, I_Cu(k), z_Cu);
    l1 = dC_dt(C_Fe2(k), C_Fe2_o, -I_Fe(k), z_Fe);
    l2 = dC_dt(C_Fe2(k)+h*l1/2, C_Fe2_o, -I_Fe(k), z_Fe);
    l3 = dC_dt(C_Fe2(k)+h*l2/2, C_Fe2_o, -I_Fe(k), z_Fe);
    l4 = dC_dt(C_Fe2(k)+h*l3, C_Fe2_o, -I_Fe(k), z_Fe);
    m1 = dC_dt(C_Fe3(k), C_Fe3_o, I_Fe(k), z_Fe);
    m2 = dC_dt(C_Fe3(k)+h*m1/2, C_Fe3_o, I_Fe(k), z_Fe);
    m3 = dC_dt(C_Fe3(k)+h*m2/2, C_Fe3_o, I_Fe(k), z_Fe);
    m4 = dC_dt(C_Fe3(k)+h*m3, C_Fe3_o, I_Fe(k), z_Fe);
    %Subplus function prevents negative concentrations in the cell
	C_Cu2(k+1) = subplus(C_Cu2(k) + 1/6*h*(k1+2*k2+2*k3+k4));
    C_Fe2(k+1) = subplus(C_Fe2(k) + 1/6*h*(l1+2*l2+2*l3+l4));
    C_Fe3(k+1) = subplus(C_Fe3(k) + 1/6*h*(m1+2*m2+2*m3+m4));
    if C_Cu2(k+1) == 0 | C_Fe2(k+1) == 0 | C_Fe3(k+1) == 0
        break
    end
	%Nernst potential calculation
	Erev_Cu(k+1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2(k+1)));
	Erev_Fe(k+1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2(k+1)/gamma_Fe3*C_Fe3(k+1));
end

eta_an = E_an - Erev_Fe;
eta_cat = E_cat - Erev_Cu;

function Y = cell_solver(I_an, E_an, E_cat, V_app, r_sol, r_hardware, Erev_Fe, Erev_Cu, S_an, S_cat)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    global temp i0_Fe alpha_Fe z_Fe i0_Cu alpha_Cu z_Cu
    eta_an = E_an - Erev_Fe;
    eta_cat = E_cat - Erev_Cu;
    Y(1) = E_an - E_cat + I_an*(r_sol+r_hardware) - V_app;
    Y(2) = I_an - i_BV(eta_an, i0_Fe, alpha_Fe, z_Fe, temp)*S_an;
    Y(3) = I_an + i_BV(eta_cat, i0_Cu, alpha_Cu, z_Cu, temp)*S_cat;
end

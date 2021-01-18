%{
	Continuous stirred tank model for base metal extraction/recovery cell
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
n = z_Cu+z_Fe;
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
%Cell Concentrations (recovery)
Ci_Cu2_cell = 0.2;
Ci_Fe2_cell = 0.5;
Ci_Fe3_cell = 0.001; 
%Bed concentrations (extraction)
Ci_Cu2_bed = 0.2;
Ci_Fe2_bed = 0.5;
Ci_Fe3_bed = 0.001;

%Electrode areas
S_cat = 100; %cm^2
S_an = 100; %cm^2
%Applied Voltage (potentiostat)
V_app = 3; %V

%Specific Surface Area of E-waste
V_bed = 0.2; %m3 Volume of bed holding the particles assuming the bed is completly full.
radius = 0.001; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
SSA = 3/radius; %m2/m3 Specific Surface area of spheres.
density = 0.6; %m3/m3 Loose packing density of equal sized spheres. Close packing density = 0.64.
S_corr = V_bed*density*SSA; %m2 Total surface area of spheres in bed.

%Mass balance on generic ion in a CSTR
%C - Concentration in CSTR, C_o - Concentration in inlet
dC_dt = @(C, C_o, I, z) ((C_o-C)*Q + I*F/z)/volume;

%Solver
tfinal = 0.001; %total time in seconds
h = 0.00001; %Step size = 1 s
t = [0:h:tfinal];

%initializing concentrations
C_Cu2_cell = zeros(1, tfinal/h+1);
C_Fe2_cell = zeros(1, tfinal/h+1);
C_Fe3_cell = zeros(1, tfinal/h+1);
C_Cu2_cell(1) = Ci_Cu2_cell;
C_Fe2_cell(1) = Ci_Fe2_cell;
C_Fe3_cell(1) = Ci_Fe3_cell;
C_Cu2_bed = zeros(1, tfinal/h+1);
C_Fe2_bed = zeros(1, tfinal/h+1);
C_Fe3_bed = zeros(1, tfinal/h+1);
C_Cu2_bed(1) = Ci_Cu2_bed;
C_Fe2_bed(1) = Ci_Fe2_bed;
C_Fe3_bed(1) = Ci_Fe3_bed;

%Nernst potential arrays
Erev_Cu_cell = zeros(1, tfinal/h+1);
Erev_Fe_cell = zeros(1, tfinal/h+1);
Erev_Cu_cell(1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2_cell(1)));
Erev_Fe_cell(1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2_cell(1)/gamma_Fe3*C_Fe3_cell(1));
Erev_Cu_bed = zeros(1, tfinal/h+1);
Erev_Fe_bed = zeros(1, tfinal/h+1);
Erev_Cu_bed(1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2_bed(1)));
Erev_Fe_bed(1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2_bed(1)/gamma_Fe3*C_Fe3_bed(1));

%Resistances
r_sol = ones(1, tfinal/h); %ohms, add calculation for this later
r_hardware = 0; %ohms

%Currents, potential, overpotential arrays
I_Cu = zeros(1, tfinal/h+1);
I_Fe = zeros(1, tfinal/h+1);
I_corr = zeros(1, tfinal/h+1);
E_an = zeros(1, tfinal/h+1);
E_cat = zeros(1, tfinal/h+1);
E_corr = zeros(1, tfinal/h+1);
%initial guess for E_corr
j0 = 0; %V

for k = 1:1:(tfinal/h)
    %solve for currents and electrode potentials in cell using fsolve
    %wrapper function for cell solver
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol(k), r_hardware, Erev_Fe_cell(k), Erev_Cu_cell(k), S_an, S_cat);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [1, 0.2, -0.2];
    x = fsolve(solver, x0);
    I_an = x(1);
    E_an(k) = x(2);
    E_cat(k) = x(3);
    I_Cu(k) = - I_an;
    I_Fe(k) = I_an;
    
    %solve for extraction/corrosion current and potential
    cor_solver = @(E_corr)cor(E_corr, [Erev_Cu_bed(k), Erev_Fe_bed(k)]);
    E_corr(k) = fzero(cor_solver, j0);
    j0 = E_corr(k);
    I_corr(k) = S_corr*i_BV(E_corr(k)-Erev_Cu_bed(k), i0_Cu, alpha_Cu, z_Cu, temp);
    
    %4th order runge-kutta to solve mass balance differential equations in
    %cell
	k1 = dC_dt(C_Cu2_cell(k), C_Cu2_bed(k), I_Cu(k), z_Cu);
	k2 = dC_dt(C_Cu2_cell(k)+h*k1/2, C_Cu2_bed(k), I_Cu(k), z_Cu);
	k3 = dC_dt(C_Cu2_cell(k)+h*k2/2, C_Cu2_bed(k), I_Cu(k), z_Cu);
	k4 = dC_dt(C_Cu2_cell(k)+h*k3, C_Cu2_bed(k), I_Cu(k), z_Cu);
    l1 = dC_dt(C_Fe2_cell(k), C_Fe2_bed(k), -I_Fe(k), z_Fe);
    l2 = dC_dt(C_Fe2_cell(k)+h*l1/2, C_Fe2_bed(k), -I_Fe(k), z_Fe);
    l3 = dC_dt(C_Fe2_cell(k)+h*l2/2, C_Fe2_bed(k), -I_Fe(k), z_Fe);
    l4 = dC_dt(C_Fe2_cell(k)+h*l3, C_Fe2_bed(k), -I_Fe(k), z_Fe);
    m1 = dC_dt(C_Fe3_cell(k), C_Fe3_bed(k), I_Fe(k), z_Fe);
    m2 = dC_dt(C_Fe3_cell(k)+h*m1/2, C_Fe3_bed(k), I_Fe(k), z_Fe);
    m3 = dC_dt(C_Fe3_cell(k)+h*m2/2, C_Fe3_bed(k), I_Fe(k), z_Fe);
    m4 = dC_dt(C_Fe3_cell(k)+h*m3, C_Fe3_bed(k), I_Fe(k), z_Fe);
    %Subplus function prevents negative concentrations in the cell
	C_Cu2_cell(k+1) = subplus(C_Cu2_cell(k) + 1/6*h*(k1+2*k2+2*k3+k4));
    C_Fe2_cell(k+1) = subplus(C_Fe2_cell(k) + 1/6*h*(l1+2*l2+2*l3+l4));
    C_Fe3_cell(k+1) = subplus(C_Fe3_cell(k) + 1/6*h*(m1+2*m2+2*m3+m4));
    
    %Solve mass balance DE on ions in extraction bed
    n1 = dC_dt(C_Cu2_bed(k), C_Cu2_cell(k), I_corr(k), z_Cu);
	n2 = dC_dt(C_Cu2_bed(k)+h*n1/2, C_Cu2_cell(k), I_corr(k), z_Cu);
	n3 = dC_dt(C_Cu2_bed(k)+h*n2/2, C_Cu2_cell(k), I_corr(k), z_Cu);
	n4 = dC_dt(C_Cu2_bed(k)+h*n3, C_Cu2_cell(k), I_corr(k), z_Cu);
    o1 = dC_dt(C_Fe2_bed(k), C_Fe2_cell(k), I_corr(k), z_Fe);
    o2 = dC_dt(C_Fe2_bed(k)+h*o1/2, C_Fe2_cell(k), I_corr(k), z_Fe);
    o3 = dC_dt(C_Fe2_bed(k)+h*o2/2, C_Fe2_cell(k), I_corr(k), z_Fe);
    o4 = dC_dt(C_Fe2_bed(k)+h*o3, C_Fe2_cell(k), I_corr(k), z_Fe);
    p1 = dC_dt(C_Fe3_bed(k), C_Fe3_cell(k), I_Fe(k), z_Fe);
    p2 = dC_dt(C_Fe3_bed(k)+h*p1/2, C_Fe3_cell(k), I_Fe(k), z_Fe);
    p3 = dC_dt(C_Fe3_bed(k)+h*p2/2, C_Fe3_cell(k), I_Fe(k), z_Fe);
    p4 = dC_dt(C_Fe3_bed(k)+h*p3, C_Fe3_cell(k), I_Fe(k), z_Fe);
    %Subplus function prevents negative concentrations in the cell
	C_Cu2_bed(k+1) = subplus(C_Cu2_bed(k) + 1/6*h*(n1+2*n2+2*n3+n4));
    C_Fe2_bed(k+1) = subplus(C_Fe2_bed(k) + 1/6*h*(o1+2*o2+2*o3+o4));
    C_Fe3_bed(k+1) = subplus(C_Fe3_bed(k) + 1/6*h*(p1+2*p2+2*p3+p4));
    
    %Check for 0 values and exit solver if any ion is depleted
    %if C_Cu2_cell(k+1) == 0 | C_Fe2_cell(k+1) == 0 | C_Fe3_cell(k+1) == 0
       % break
    %end
    %if C_Cu2_bed(k+1) == 0 | C_Fe2_bed(k+1) == 0 | C_Fe3_bed(k+1) == 0
       % break
    %end
    
	%Nernst potential calculations
	Erev_Cu_cell(k+1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2_cell(k+1)));
	Erev_Fe_cell(k+1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2_cell(k+1)/gamma_Fe3*C_Fe3_cell(k+1));
    Erev_Cu_bed(k+1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2_bed(k+1)));
    Erev_Fe_bed(k+1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2_bed(k+1)/gamma_Fe3*C_Fe3_bed(k+1));
end

eta_an = E_an - Erev_Fe_cell;
eta_cat = E_cat - Erev_Cu_cell;

function Y = cell_solver(I_an, E_an, E_cat, V_app, r_sol, r_hardware, Erev_Fe, Erev_Cu, S_an, S_cat)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    global temp i0_Fe alpha_Fe z_Fe i0_Cu alpha_Cu z_Cu
    eta_an = E_an - Erev_Fe;
    eta_cat = E_cat - Erev_Cu;
    Y(1) = E_an - E_cat + I_an*(r_sol+r_hardware) - V_app;
    Y(2) = I_an - i_BV(eta_an, i0_Fe, alpha_Fe, z_Fe, temp)*S_an;
    Y(3) = I_an + i_BV(eta_cat, i0_Cu, alpha_Cu, z_Cu, temp)*S_cat;
end

function [func] = cor(Ecorr, Erev)
  %1 = Cu, 2= Fe
  global i0_Cu i0_Fe alpha_Cu alpha_Fe z_Cu z_Fe temp
  i_Cu = i_BV(Ecorr-Erev(1), i0_Cu, alpha_Cu, z_Cu, temp);
  i_Fe = i_BV(Ecorr-Erev(2), i0_Fe, alpha_Fe, z_Fe, temp);
  func = i_Cu + i_Fe;
  end


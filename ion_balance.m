function dt = ion_balance(t, Cm, temp, pres, vol_cell, vol_lch, Q, S_an, S_cat, V_app, n_particles, l, A_cell)
    %t = time span (s)
    %Cm = ionic concentration and solid mass vector (M, kg)
    %Concentration order: Cu2 Sn2 Al3 Pb2 Fe2 Fe3 Ag Au Pd2 H Cl
    %Solid mass order: Inert Cu Sn Al Pb Fe Ag Au Pd
    %pres = system pressure (atm)
    %vol_cell = electrowinning cell volume (L)
    %vol_lch = leaching vessel volume (L)
    %Q = circulation flowrate (L/s)
    %S_an = cell anodic surface area (cm^2);
    %Sn_cat = cell cathodic surface area (cm^2);
    %V_app = applied voltage to cell (V);
    %n_particles = number of PCB particles as calculated by radius/volume;
    %l = cell distance between electrodes;
    %A_cell = cell cross sectional area;
    
    propertiesBaseMetals
    
    %calculate total mass and wt fractions based on partial masses at
    %timestep
    m_PCB = (Cm(23:31).');
    m_PCB_total = sum(m_PCB);
    wtfrac_PCB = m_PCB/m_PCB_total;
    %Convert weight to volume and volume fraction
    V_PCB = m_PCB./rho; %Element volumes in m^3
    vfrac_PCB = V_PCB/sum(V_PCB);
    V_PCB_total = sum(V_PCB);
    
    %calculate particle radius in m
    r_particles = (V_PCB_total/(n_particles*4*pi)*3)^(1/3);
    
    %Surface area calculation for corrosion
    packing_density = 0.6;
    SSA = 3/r_particles; %m2/m3 Specific Surface area of spheres.
    S_PCB_total = V_PCB_total*packing_density*SSA; %m2 Total surface area of spheres in lch.
    S_PCB = vfrac_PCB*S_PCB_total; %array with exposed surface area of each component
    S_PCB = subplus(S_PCB*100^2); %convert m^2 to cm^2
    
    [Erev_cell, Erev_lch] = nernstPotential(Cm,temp,aH2,aO2);
    
    %resistance calculation for IR drop in cell
    kappa = 1000*(Cm(1)*lamda(1)+Cm(2)*lamda(2)+Cm(3)*lamda(3)+...
        Cm(4)*lamda(4)+Cm(5)*lamda(5)+Cm(6)*lamda(6)+Cm(7)*lamda(7)+...
        Cm(8)*lamda(8)+Cm(9)*lamda(9)+Cm(10)*lamda(10)+Cm(11)*lamda(11));
    r_sol = l/A_cell/kappa*100;
    r_hardware = 1; %ohms
    
    %%%Electrowinning Cell solving%%%
    %solve cell currents and electrode potentials
    onCathode = [1 1 1 1 0 1 1 1 1 1 0];
    onAnode = -(onCathode-1);
    solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol, r_hardware,...
        Erev_cell, onCathode, S_an, S_cat, temp);
    %initial guesses [I_an, E_an, E_cat]
    x0 = [0.2, 0.1, -0.1];
    options = optimoptions(@fsolve, 'Display','final', 'MaxFunctionEvaluations', 3000);
    x = fsolve(solver, x0, options);
    I = x(1);
    E_an = x(2);
    E_cat = x(3);
    eta_cat = E_cat - Erev_cell;
    eta_an = E_an - Erev_cell;
    i_cat = onCathode.*i_BV(eta_cat, i0, alphas, z, temp);
    I_cat = i_cat*S_cat;
    i_an = onAnode.*i_BV(eta_an, i0, alphas, z, temp);
    I_an = i_an*S_an;
    I_cell = I_cat+I_an; %overall current for rxn i in cell
    
    %%%Leaching Unit solving%%%
    %solve extraction lch corrosion rate
    j0 = Erev_lch(5)*0.90; %Initial guess for E_corr, V
    cor_solver = @(E_corr)cor(E_corr, Erev_lch, S_PCB, temp);
    E_corr = fzero(cor_solver, j0);
    i_corr = i_BV(E_corr-Erev_lch, i0, alphas, z, temp);
    S_corr = [S_PCB(2:5) sum(S_PCB) S_PCB(6:9) 0 sum(S_PCB)];
    I_corr = S_corr.*i_corr;
    
    %%%Solving for mass balance differentials%%%%
    %Calculate concentration/mass balances
    dt = zeros(size(Cm));
    %Electrowinning cell concentration balances
    dt(1) = ((Cm(12)-Cm(1))*Q + I_cell(1)/F/z(1))/vol_cell; %Cu2+
    dt(2) = ((Cm(13)-Cm(2))*Q + I_cell(2)/F/z(2))/vol_cell; %Sn2+
    dt(3) = ((Cm(14)-Cm(3))*Q + I_cell(3)/F/z(3))/vol_cell; %Al3+
    dt(4) = ((Cm(15)-Cm(4))*Q + I_cell(4)/F/z(4))/vol_cell; %Pb2+
    dt(5) = ((Cm(16)-Cm(5))*Q + I_cell(6)/F/z(6)-I_cell(5)/F/z(5))/vol_cell; %Fe2+
    dt(6) = ((Cm(17)-Cm(6))*Q + I_cell(5)/F/z(6))/vol_cell; %Fe3+
    dt(7) = ((Cm(18)-Cm(7))*Q + I_cell(7)/F/z(7))/vol_cell; %Ag+
    dt(8) = ((Cm(19)-Cm(8))*Q + I_cell(8)/F/z(8))/vol_cell; %Au+
    dt(9) = ((Cm(20)-Cm(9))*Q + I_cell(9)/F/z(9))/vol_cell; %Pd2+
    dt(10) = ((Cm(21)-Cm(10))*Q + I_cell(10)/F/z(10) + I_cell(11)/F/z(11))/vol_cell; %H+
    dt(11) = (Cm(22)-Cm(11))*Q/vol_cell; %Cl-
    
    %Leaching vessel concentration balances
    dt(12) = ((Cm(1)-Cm(12))*Q + I_corr(1)/F/z(1))/vol_lch; %Cu2+
    dt(13) = ((Cm(2)-Cm(13))*Q + I_corr(2)/F/z(2))/vol_lch; %Sn2+
    dt(14) = ((Cm(3)-Cm(14))*Q + I_corr(3)/F/z(3))/vol_lch; %Al3+
    dt(15) = ((Cm(4)-Cm(15))*Q + I_corr(4)/F/z(4))/vol_lch; %Pb2+
    dt(16) = ((Cm(5)-Cm(16))*Q + I_corr(6)/F/z(6)-I_corr(5)/F/z(5))/vol_lch; %Fe2+
    dt(17) = ((Cm(6)-Cm(17))*Q + I_corr(5)/F/z(6))/vol_lch; %Fe3+
    dt(18) = ((Cm(7)-Cm(18))*Q + I_corr(7)/F/z(7))/vol_lch; %Ag+
    dt(19) = ((Cm(8)-Cm(19))*Q + I_corr(8)/F/z(8))/vol_lch; %Au+
    dt(20) = ((Cm(9)-Cm(20))*Q + I_corr(9)/F/z(9))/vol_lch; %Pd2+
    dt(21) = ((Cm(10)-Cm(21))*Q + I_corr(10)/F/z(10) + I_corr(11)/F/z(11))/vol_lch; %H+
    dt(22) = (Cm(11)-Cm(22))*Q/vol_lch; %Cl-
    
    %PCB metal mass balances
    dt(23) = 0; %Inert material 
    dt(24:31) = -mw(2:9).*[I_corr(1:4) I_corr(6:9)]/F./[z(1:4) z(6:9)]; %Solid metals Cu to Pd
    %display current time step
    t
end

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
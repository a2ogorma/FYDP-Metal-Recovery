function results = BaseMetalCell(initSet,paramSet)
    propertiesBaseMetals %calls all the property values as set in the applicable file. Note this is currently only applicable to base metals
    %{
	Continuous stirred tank model for base metal extraction/recovery cell
    Using ode45
	To Add: Applied Voltage calc. Non-constant volume?
    Order of compounds: Inert Cu Sn Al Pb Fe Zn Ca Ni Ag Au Pd
    %}
    tic
    % system parameters
    temp = paramSet.temp; %K
    pres = paramSet.pres; % atm
    vol_cell = paramSet.vol_cell; %L
    Q = paramSet.Q;%; % L/s (flowrate)
    %Electrode areas
    S_cat = paramSet.S_cat; %cm^2
    S_an = paramSet.S_an; %cm^2
    %Cross sectional area of cell
    A_cell = paramSet.A_cell; %cm^2
    %Length b/w electrodes
    l = paramSet.l; %cm
    %Applied Voltage (potentiostat)
    V_app = paramSet.V_app; %V
    %Extraction vessel parameters
    vol_lch = paramSet.vol_lch; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
    tfinal = paramSet.tfinal; %s

    %PCB characteristics
    m_PCB_total = initSet.solidPCB.m_PCB_total; %kg Mass of crushed PCBs
    r_particles = initSet.solidPCB.r_particles; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
    %Weight fraction composition of PCB
    %Inert Cu Sn Al Pb Fe 
    wtfrac_PCB = initSet.solidPCB.wtfrac_PCB;
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

    %initializing solution concentration and solid mass vector
    Cm_i = [initSet.solution.Ci_cell initSet.solution.Ci_lch m_PCB];

    %solve concentration profiles
    tspan = [0 tfinal];
    options = odeset('NonNegative',1:31);
    balance_solver = @(t, Cm) ion_balance(t, Cm, temp, pres, vol_cell, vol_lch, Q, S_an, S_cat, V_app, n_particles, l, A_cell);
    [t, Cm] = ode15s(balance_solver, tspan, Cm_i, options);

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
        onCathode = [1 1 1 1 1 1 1 1 1 1 0];
        onAnode = [0 0 0 0 1 0 0 0 0 0 1];
        solver = @(x) cell_solver(x(1), x(2), x(3), V_app, r_sol(j), r_hardware, ...
            Erev_cell(j,:), onCathode, onAnode, S_an, S_cat, temp);
        %initial guesses [I_an, E_an, E_cat]
        x0 = [0.2, 0.1, -0.1];
        options = optimoptions(@fsolve, 'Display','final', 'MaxFunctionEvaluations', 3000);
        x = fsolve(solver, x0, options);
        I(j) = x(1);
        E_an(j) = x(2);
        E_cat(j) = x(3);
        eta_cat(j,:) = E_cat(j) - Erev_cell(j,:);
        eta_an(j,:) = E_an(j) - Erev_cell(j,:);
        i_cat(j,:) = onCathode.*i_BV(eta_cat(j,:), i0, iL, alphas, z, temp);
        I_cat(j,:) = i_cat(j,:)*S_cat;
        i_an(j,:) = onAnode.*i_BV(eta_an(j,:), i0, iL, alphas, z, temp);
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
        on_PCB_cathode = [1 1 1 1 1 1 1 1 1 1 1];
        on_PCB_anode = [1 1 1 1 1 1 1 1 1 0 1];
        j0 = 0; %Initial guess for E_corr, V
        cor_solver = @(E_corr)cor(E_corr, Erev_lch(j,:), S_PCB(j,:), ...
            on_PCB_cathode, on_PCB_anode, temp);
        E_corr(j) = fzero(cor_solver, j0);
        i_BV_corr = i_BV(E_corr(j)-Erev_lch(j,:), i0, iL, alphas, z, temp);
        i_corr_an = on_PCB_anode.*subplus(i_BV_corr);
        i_corr_cat = on_PCB_cathode.*(-subplus(-i_BV_corr));
        i_corr(j,:) = i_corr_an + i_corr_cat;
        S_corr(j,:) = [S_PCB(j,2:5) sum(S_PCB(j,2:9)) S_PCB(j,6:9) sum(S_PCB(j,2:9)) sum(S_PCB(j,2:9))];
        I_corr(j,:) = S_corr(j,:).*i_corr(j,:);
    end
    
    %% Setting up results output struct
    results.init.paramSet = paramSet;
    results.init.initSet = initSet;
    results.init.n_particles = n_particles;
    
    results.leaching.solutionConc = Cm(:,12:22);
    results.leaching.E_corr = E_corr;
    results.leaching.Erev_lch = Erev_lch;
    results.leaching.i_corr = i_corr;
    results.leaching.I_corr = I_corr;
    
    results.electrowinning.solutionConc = Cm(:,1:11);
    results.electrowinning.E_an = E_an;
    results.electrowinning.E_cat = E_cat;
    results.electrowinning.Erev_cell = Erev_cell;
    results.electrowinning.eta_an = eta_an;
    results.electrowinning.eta_cat = eta_cat;
    results.electrowinning.I = I;
    results.electrowinning.I_an = I_an;
    results.electrowinning.i_an = i_an;
    results.electrowinning.I_cat = I_cat;
    results.electrowinning.i_cat = i_cat;
    results.electrowinning.I_cell = I_cell;
    results.electrowinning.r_sol = r_sol;
    
    results.PCB.r_particles = r_particles;
    results.PCB.S_corr = S_corr;
    results.PCB.S_PCB = S_PCB;
    results.PCB.S_PCB_total = S_PCB_total;
    results.PCB.SSA = SSA;
    results.PCB.massRem = m_PCB;
    results.PCB.massTotal = m_PCB_total;
    results.PCB.V_PCB = V_PCB;
    results.PCB.V_PCB_total = V_PCB_total;
    results.PCB.vfrac_PCB = vfrac_PCB;
    results.PCB.wtfrac_PCB = wtfrac_PCB;
    massStart = m_PCB(1,:);
    results.PCB.pctRem = results.PCB.massRem./massStart;
    results.PCB.massRec.Copper(1) = 0;
    results.PCB.massRec.Tin(1) = 0;
    results.PCB.massRec.Aluminum(1) = 0;
    results.PCB.massRec.Lead(1) = 0;
    results.PCB.massRec.Iron(1) = 0;
    results.PCB.massRec.Silver(1) = 0;
    results.PCB.massRec.Gold(1) = 0;
    results.PCB.massRec.Palladium(1) = 0;
    for j = 2:1:length(I_cell)
        I_cell(1:j,1)
        results.PCB.massRec.Copper(j) = -mw(2)*trapz(t(1:j), I_cat(1:j,1))/F/z(1)/1000;
        results.PCB.massRec.Tin(j) = -mw(3)*trapz(t(1:j), I_cat(1:j,2))/F/z(2)/1000;
        results.PCB.massRec.Aluminum(j) = -mw(4)*trapz(t(1:j), I_cat(1:j,3))/F/z(3)/1000;
        results.PCB.massRec.Lead(j) = -mw(5)*trapz(t(1:j), I_cat(1:j,4))/F/z(4)/1000;
        results.PCB.massRec.Iron(j) = -mw(6)*trapz(t(1:j), I_cat(1:j,6))/F/z(6)/1000;
        results.PCB.massRec.Silver(j) = -mw(7)*trapz(t(1:j), I_cat(1:j,7))/F/z(7)/1000;
        results.PCB.massRec.Gold(j) = -mw(8)*trapz(t(1:j), I_cat(1:j,8))/F/z(8)/1000;
        results.PCB.massRec.Palladium(j) = -mw(9)*trapz(t(1:j), I_cat(1:j,9))/F/z(9)/1000;
    end
    results.PCB.pctRec.Copper = -results.PCB.massRec.Copper./massStart(2);
    results.PCB.pctRec.Tin = -results.PCB.massRec.Tin./massStart(3);
    results.PCB.pctRec.Aluminum = -results.PCB.massRec.Aluminum./massStart(4);
    results.PCB.pctRec.Lead = -results.PCB.massRec.Lead./massStart(5);
    results.PCB.pctRec.Iron = -results.PCB.massRec.Iron./massStart(6);
    results.PCB.pctRec.Silver = -results.PCB.massRec.Silver./massStart(7);
    results.PCB.pctRec.Gold = -results.PCB.massRec.Gold./massStart(8);
    results.PCB.pctRec.Palladium = -results.PCB.massRec.Palladium./massStart(9);
    
    results.t = t;
    results.simulationTime = toc;
end

%%corrosion function
function func = cor(Ecorr, Erev, S_PCB, on_PCB_cat, on_PCB_an, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 iL alphas z
    i_BV_corr = i_BV(Ecorr-Erev, i0, iL, alphas, z, temp);
    i_corr_an = on_PCB_an.*subplus(i_BV_corr);
    i_corr_cat = on_PCB_cat.*(-subplus(-i_BV_corr));
    i_corr = i_corr_an + i_corr_cat;
    %arrange surface areas in proper order
    S_corr = [S_PCB(2:5) sum(S_PCB(2:9)) S_PCB(6:9) sum(S_PCB(2:9)) sum(S_PCB(2:9))];
    if abs(sum(S_corr.*i_corr)) == Inf
        disp(S_corr)
        disp(i_corr)
        error("Corrosion currents infinite");
    end
    func = sum(S_corr.*i_corr);
end
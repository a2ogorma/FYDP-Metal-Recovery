function results = BaseMetalCell(initSet,paramSet)
    %{
	Continuous stirred tank model for base metal extraction/recovery cell
    Order of compounds: Inert Cu Sn Al Pb Fe Ag Au Pd
    %}
    tic
    
    %% Solving concentration/mass matrix Cm using ODE solvers
    m_deposited = initSet.m_deposited;
    %for solution, 1 is for Cl- base metal system, 2 is for S2o3- precious metal system
    solution = initSet.solution.type;
    propertiesMetals %calls all the property values as set in the applicable file.
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
    %Operation mode 
    mode = paramSet.mode; % 1 - potentiostat, 2 - galvanostat
    %Applied Voltage (potentiostat)
    V_app = paramSet.V_app; %V
    %Applied Current (galvanostat)
    I_app = paramSet.I_app; %A
    %Extraction vessel parameters
    vol_lch = paramSet.vol_lch; %L (Initial) volume of bed holding the particles assuming the bed is completly full.
    tfinal = paramSet.tfinal; %s
    
    %fsolve options
    foptions = paramSet.foptions;
    
    %PCB characteristics
    m_PCB_total = initSet.solidPCB.m_PCB_total; %kg Mass of crushed PCBs
    r_particles = initSet.solidPCB.r_particles; %m Radius of particles. Must be 2.873 (or greater) times smaller than the radius of the cylinder.
    %Weight fraction composition of PCB
    %Inert Cu Sn Fe Ag Au Pd
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
    S_PCB_total = V_PCB_total*SSA; %m2 Total surface area of spheres in lch.
    S_PCB = vfrac_PCB*S_PCB_total; %array with exposed surface area of each component
    S_PCB = S_PCB*100^2; %convert m^2 to cm^2

    %Calculate number of crushed PCB particles
    n_particles = sum(V_PCB)*3/(4*pi*r_particles^3);

    %initializing solution concentration and solid mass vector
    Cm_i = [initSet.solution.Ci_cell initSet.solution.Ci_cell initSet.solution.Ci_lch m_PCB m_deposited];

    %Function for detect discontinuities in corrosion and solver failures
    if mode == 1
        discont_event = @(t,Cm) discont(t, Cm, temp, pres, vol_cell, ...
            vol_lch, Q, S_an, S_cat, mode, V_app, n_particles, l, A_cell, ...
            solution, foptions);
    elseif mode == 2
        discont_event = @(t,Cm) discont(t, Cm, temp, pres, vol_cell, ... 
            vol_lch, Q, S_an, S_cat, mode, I_app, n_particles, l, A_cell, ...
            solution, foptions);
    end
    %}
    %solve concentration profiles
    tspan = [0 tfinal];
    %options = odeset('NonNegative', 1:41);
    %options = odeset('Events', discont_event);
    options = odeset('NonNegative',1:41, 'Events', discont_event);
    if mode == 1
        balance_solver = @(t, Cm) ion_balance(t, Cm, temp, pres, vol_cell, ...
            vol_lch, Q, S_an, S_cat, mode, V_app, n_particles, l, A_cell, ...
            solution, foptions);
    elseif mode == 2
        balance_solver = @(t, Cm) ion_balance(t, Cm, temp, pres, vol_cell, ...
            vol_lch, Q, S_an, S_cat, mode, I_app, n_particles, l, A_cell, ...
            solution, foptions);
    end
    [t, Cm, te, Cme, ie] = ode15s(balance_solver, tspan, Cm_i, options);
    %[t, Cm] = ode45(balance_solver, tspan, Cm_i,options);
    
    disp("Post-processing")
    
    %% Back calculating currents and potentials using solution for Cm matrix:
    for j = 1:1:length(t)
        CmStep = Cm(j,:);
        %Nernst Potentials
        [Erev_cat(j,:), Erev_an(j,:), Erev_lch(j,:)] = nernstPotential(CmStep,temp,solution);
        %resistance calculation for IR drops
        kappa_cat = 1000*(CmStep(1)*lamda(1)+CmStep(2)*lamda(2)+CmStep(3)*lamda(3)+...
            CmStep(4)*lamda(4)+CmStep(5)*lamda(5)+CmStep(6)*lamda(6)+CmStep(7)*lamda(7)+...
            CmStep(8)*lamda(8)+CmStep(9)*lamda(9)+CmStep(10)*lamda(10));
        kappa_an = 1000*(CmStep(11)*lamda(1)+CmStep(12)*lamda(2)+CmStep(13)*lamda(3)+...
            CmStep(14)*lamda(4)+CmStep(15)*lamda(5)+CmStep(16)*lamda(6)+CmStep(17)*lamda(7)+...
            CmStep(18)*lamda(8)+CmStep(19)*lamda(9)+CmStep(20)*lamda(10));
        kappa_avg = mean([kappa_an kappa_cat]);
        r_sol = l/A_cell/kappa_avg*100;
        r_hardware = 10; %ohms
        
        %cell volume division
        vol_cat = vol_cell/2;
        vol_an = vol_cat;
        
        %surface area calculation for cathode
        v_cat = CmStep(38:43)./rho(2:7);
        vfrac_cat = v_cat/sum(v_cat);
        S_cat_p(j,:) = (S_cat*vfrac_cat);
        
        %limiting current calculations for cathode side
        iLc_cat(j,1) = z(1)*F*km(1)*CmStep(1)+eps;
        iLc_cat(j,2) = z(2)*F*km(2)*CmStep(2)+eps;
        iLc_cat(j,3) = z(3)*F*km(4)*CmStep(4)+eps;
        iLc_cat(j,4) = z(4)*F*km(3)*CmStep(3)+eps;
        iLc_cat(j,6) = -1;
        iLc_cat(j,8) = z(8)*F*km(10)*CmStep(10)+eps;
        iLc_cat(j,10) = z(10)*F*km(8)*CmStep(8)+eps;
        iLc_cat(j,11) = z(11)*F*km(8)*CmStep(8)+eps;

        iLa_cat(j,:) = -1*ones(1,11);
        iLa_cat(j,3) = z(5)*F*km(3)*CmStep(3)+eps;
        iLa_cat(j,6) = z(6)*F*km(9)*CmStep(9)+eps;
        iLa_cat(j,8) = z(8)*F*km(9)*CmStep(9)/4+eps;

        if solution == 1 %Cl-, base metal system
            iLc_cat(j,5) = z(5)*F*km(5)*CmStep(5)+eps;
            iLc_cat(j,7) = z(7)*F*km(6)*CmStep(6)+eps;
            iLc_cat(j,9) = z(9)*F*km(7)*CmStep(7)+eps;
        else %S2O3, precious metal system
            iLc_cat(j,5) = z(5)*F*km(5)*CmStep(5)+eps;
            iLc_cat(j,7) = z(7)*F*km(6)*CmStep(6)+eps;
            iLc_cat(j,9) = z(9)*F*km(7)*CmStep(7)+eps;
            iLa_cat(j,5) = z(5)*F*km(9)*CmStep(9)/2+eps;
            iLa_cat(j,7) = z(7)*F*km(9)*CmStep(9)/2+eps;
            iLa_cat(j,9) = z(9)*F*km(9)*CmStep(9)/4+eps;
        end

        %limiting current calculations for anode side
        iLc_an(j,1) = z(1)*F*km(1)*CmStep(11)+eps;
        iLc_an(j,2) = z(2)*F*km(2)*CmStep(12)+eps;
        iLc_an(j,3) = z(3)*F*km(4)*CmStep(14)+eps;
        iLc_an(j,4) = z(4)*F*km(3)*CmStep(13)+eps;
        iLc_an(j,6) = -1;
        iLc_an(j,8) = z(8)*F*km(10)*CmStep(20)+eps;
        iLc_an(j,10) = z(10)*F*km(8)*CmStep(18)+eps;
        iLc_an(j,11) = z(11)*F*km(8)*CmStep(18)+eps;

        iLa_an(j,:) = -1*ones(1,11);
        iLa_an(j,3) = z(5)*F*km(3)*CmStep(13)+eps;
        iLa_an(j,6) = z(6)*F*km(9)*CmStep(19)+eps;
        iLa_an(j,8) = z(8)*F*km(9)*CmStep(19)/4+eps;

        if solution == 1 %Cl-, base metal system
            iLc_an(j,5) = z(5)*F*km(5)*CmStep(15)+eps;
            iLc_an(j,7) = z(7)*F*km(6)*CmStep(16)+eps;
            iLc_an(j,9) = z(9)*F*km(7)*CmStep(17)+eps;
        else %S2O3, precious metal system
            iLc_an(j,5) = z(5)*F*km(5)*CmStep(15)+eps;
            iLc_an(j,7) = z(7)*F*km(6)*CmStep(16)+eps;
            iLc_an(j,9) = z(9)*F*km(7)*CmStep(17)+eps;
            iLa_an(j,5) = z(5)*F*km(9)*CmStep(19)/2+eps;
            iLa_an(j,7) = z(7)*F*km(9)*CmStep(19)/2+eps;
            iLa_an(j,9) = z(9)*F*km(9)*CmStep(19)/4+eps;
        end
        
        foptions = optimoptions(@fsolve, 'Display','off', ...
        'MaxFunctionEvaluations', 4000);
        %%%Electrowinning Cell solving%%%
        %solve cell currents and electrode potentials
        onCathode = [1 1 1 1 1 0 1 1 1 1 1];
        onAnode = [0 0 1 0 0 0 0 0 0 0 1];
        if mode == 1
            solver = @(x) cell_solver_p(x(1), x(2), x(3), V_app, r_sol, r_hardware, ...
                Erev_cat(j,:), Erev_an(j,:), iLa_cat(j,:), iLc_cat(j,:), iLa_an(j,:), iLa_cat(j,:), onCathode, ...
                onAnode, S_an, S_cat_p(j,:), temp);
            %initial guesses [I_an, E_an, E_cat]
            x0 = [0.2, 0.1, -0.1];
            x = fsolve(solver, x0, foptions);
            I = x(1);
        elseif mode == 2
            solver = @(x) cell_solver_g(x(1), x(2), x(3), I_app, r_sol, r_hardware, ...
                Erev_cat(j,:), Erev_an(j,:), iLa_cat(j,:), iLc_cat(j,:), iLa_an(j,:), iLa_cat(j,:), onCathode, ...
                onAnode, S_an, S_cat_p(j,:), temp);
            %initial guesses [V, E_an, E_cat]
            x0 = [2, 0.5, -0.5];
            x = fsolve(solver, x0, foptions);
            V = x(1);
        end
        E_an(j) = x(2);
        E_cat(j) = x(3);
        eta_cat(j,:) = E_cat(j) - Erev_cat(j,:);
        eta_an(j,:) = E_an(j) - Erev_an(j,:);
        i_cat(j,:) = onCathode.*(i_BV(eta_cat(j,:), i0, iLa_cat(j,:), iLc_cat(j,:), alphas, z, temp));
        i_cat_cat = -subplus(-i_cat(j,:));
        i_cat_an = subplus(i_cat(j,:));
        I_cat_cat = i_cat_cat*S_cat;
        I_cat_cat(6) = 0; %silver can't be plated from AgCl(s)
        I_cat_an = i_cat_an.*[S_cat_p(j,1:2) S_cat S_cat_p(j,3:4) S_cat_p(j,4:5) S_cat_p(j,5:6) 0 S_cat];
        I_cat(j,:) = I_cat_cat+I_cat_an;
        i_an(j,:) = onAnode.*i_BV(eta_an(j,:), i0, iLa_an(j,:), iLc_an(j,:), alphas, z, temp);
        I_an(j,:) = i_an(j,:)*S_an;
        I_cell(j,:) = I_cat(j,:)+I_an(j,:); %overall current for rxn i in cell
        I_cell_err(j) = sum(I_cell(j,:));
        
        m_PCB(j,:) = (CmStep(31:37));
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
        S_PCB_total(j) = V_PCB_total(j)*SSA(j); %m2 Total surface area of spheres in lch.
        S_PCB(j,:) = vfrac_PCB(j,:)*S_PCB_total(j); %array with exposed surface area of each component
        S_PCB(j,:) = S_PCB(j,:)*100^2; %convert m^2 to cm^2

        %%%Leaching Unit solving%%%
        iLc_corr(j,1) = z(1)*F*km(1)*CmStep(21)+eps;
        iLc_corr(j,2) = z(2)*F*km(2)*CmStep(22)+eps;
        iLc_corr(j,3) = z(3)*F*km(4)*CmStep(24)+eps;
        iLc_corr(j,4) = z(4)*F*km(3)*CmStep(23)+eps;
        iLc_corr(j,6) = -1;
        iLc_corr(j,8) = z(8)*F*km(10)*CmStep(30)+eps;
        iLc_corr(j,10) = z(10)*F*km(8)*CmStep(18)+eps;
        iLc_corr(j,11) = z(11)*F*km(8)*CmStep(18)+eps;

        iLa_corr(j,:) = -1*ones(1,11);
        iLa_corr(j,3) = z(3)*F*km(3)*CmStep(23)+eps;

        if solution == 1
            iLc_corr(j,5) = z(5)*F*km(5)*CmStep(25)+eps;
            iLc_corr(j,7) = z(7)*F*km(6)*CmStep(26)+eps;
            iLc_corr(j,9) = z(9)*F*km(7)*CmStep(27)+eps;
        else
            iLc_corr(j,5) = z(5)*F*km(5)*CmStep(25)+eps;
            iLc_corr(j,7) = z(7)*F*km(6)*CmStep(26)+eps;
            iLc_corr(j,9) = z(9)*F*km(7)*CmStep(27)+eps;
            iLa_corr(j,5) = z(5)*F*km(9)*CmStep(29)/2+eps;
            iLa_corr(j,7) = z(7)*F*km(9)*CmStep(29)/2+eps;
            iLa_corr(j,9) = z(9)*F*km(9)*CmStep(29)/4+eps;
        end
        
        %solve extraction lch corrosion rate
        on_PCB_cathode = [1 1 1 1 1 0 1 1 1 1 1];
        on_PCB_anode = [1 1 1 1 1 1 1 1 1 0 1];
        j0 = (Erev_lch(j,1)+Erev_lch(j,5))/2; %Initial guess for E_corr, V
        cor_solver = @(E_corr)cor(E_corr, Erev_lch(j,:), iLa_corr(j,:), iLc_corr(j,:), S_PCB, ... 
            on_PCB_cathode, on_PCB_anode, temp);
        %E_corr = fminbnd(cor_solver,min(Erev_lch),max(Erev_lch));
        E_corr(j) = fsolve(cor_solver, j0, foptions);
        i_BV_corr = i_BV(E_corr(j)-Erev_lch(j,:), i0, iLa_corr(j,:), iLc_corr(j,:), alphas, z, temp);
        i_corr_an = on_PCB_anode.*subplus(i_BV_corr);
        i_corr_cat = on_PCB_cathode.*(-subplus(-i_BV_corr));
        i_corr(j,:) = i_corr_an + i_corr_cat;
        %arrange surface areas in proper order
        S_corr(j,:) = [S_PCB(2:3) sum(S_PCB(2:7)) S_PCB(4:5) S_PCB(5:6) S_PCB(6:7) sum(S_PCB(2:7)) sum(S_PCB(2:7))];
        I_corr(j,:) = S_corr(j,:).*i_corr(j,:);
        I_corr_err(j) = sum(I_corr(j,:));
        
        %dCm_dt(j,:) = balance_solver(t(j), CmStep');
    end
    
    %% Setting up results output struct
    results.init.paramSet = paramSet;
    results.init.initSet = initSet;
    results.init.n_particles = n_particles;
    results.Cm = Cm;
    results.t = t;
    results.te = te;
    results.Cme = Cme;
    results.ie = ie;
    %results.dCm_dt = dCm_dt;
    
    results.leaching.C = Cm(:,21:30);
    results.leaching.E_corr = E_corr;
    results.leaching.Erev_lch = Erev_lch;
    results.leaching.i_corr = i_corr;
    results.leaching.I_corr = I_corr;
    results.leaching.I_corr_err = I_corr_err;
    results.leaching.iLa_corr = iLa_corr;
    results.leaching.iLc_corr = iLc_corr;
    
    results.electrowinning.catholyte_C = Cm(:,1:10);
    results.electrowinning.anolyte_C = Cm(:,11:20);
    results.electrowinning.E_an = E_an;
    results.electrowinning.E_cat = E_cat;
    results.electrowinning.Erev_cat = Erev_cat;
    results.electrowinning.Erev_an = Erev_an;
    results.electrowinning.eta_an = eta_an;
    results.electrowinning.eta_cat = eta_cat;
    results.electrowinning.I_an = I_an;
    results.electrowinning.i_an = i_an;
    results.electrowinning.I_cat = I_cat;
    results.electrowinning.i_cat = i_cat;
    results.electrowinning.I_cell = I_cell;
    results.electrowinning.I_cell_err = I_cell_err;
    results.electrowinning.iLa_cat = iLa_cat;
    results.electrowinning.iLc_cat = iLc_cat;
    results.electrowinning.iLa_an = iLa_an;
    results.electrowinning.iLc_an = iLc_an;
    results.electrowinning.r_sol = r_sol;
    results.electrowinning.S_cat_p = S_cat_p;
    results.electrowinning.m_plated = Cm(:,38:43);
    
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
    results.simulationTime = toc;
end

%{
%% corrosion function
function func = cor(Ecorr, Erev, iL, S_PCB, on_PCB_cat, on_PCB_an, psbl_lch, temp)
    %units: I [A], E [V], V_app [V], r [ohms], S [cm^2]
    %Erev: Array of nernst potentials
    global i0 alphas z
    i_BV_corr = i_BV(Ecorr-Erev, i0, iL, alphas, z, temp);
    i_corr_an = on_PCB_an.*subplus(i_BV_corr);
    i_corr_cat = psbl_lch.*on_PCB_cat.*(-subplus(-i_BV_corr));
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
%}
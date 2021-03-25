clear all
sims = 10;
for run = 1:1:sims
    %% Initial condition specifications
    solution = 1; %1 is Cl- base metal, 2 is S2O3 precious metal
    propertiesMetals;
    initSetBase = struct;
    %characteristics of solid PCB input
    initSetBase.solidPCB.wtfrac_PCB = [0.0845266104991119	0.815039475438435	0.0971570235621975	0.00210776765005768	0.000701689614615871	0.000313061520367081	0.000154371715215492];
    %assuming cycle time is 20 hrs
    tfinal = 50*3600;
    %Assuming 100000 kg/yr waste input
    initSetBase.solidPCB.m_PCB_total = 100000/(8760*0.91)*tfinal/3600;
    global rho
    V_PCB_total = sum(initSetBase.solidPCB.m_PCB_total.*initSetBase.solidPCB.wtfrac_PCB./rho)*1000;%L
    %Particle radius, m
    initSetBase.solidPCB.r_particles = randintF(1,10,0)/1000; 

    %characteristics of starting solution
    initSetBase.solution.type = solution;%1 is Cl- base metal, 2 is S2O3 precious metal
    %initial concentrations in mol/L
    %Cell Concentrations (recovery)
    initSetBase.solution.Ci_Cu2_cell = 0.0001;
    initSetBase.solution.Ci_Sn2_cell = 0.0001;
    initSetBase.solution.Ci_Fe2_cell = 0.001;
    initSetBase.solution.Ci_Fe3_cell = randintF(0.4,1,0);
    initSetBase.solution.Ci_Ag_cell = 0.00;
    initSetBase.solution.Ci_Au3_cell = 0.0;
    initSetBase.solution.Ci_Pd2_cell = 0.0;
    initSetBase.solution.Ci_H_cell = 0.5;
    %Calculation to ensure electrolyte has net neutral charge
    initSetBase.solution.Ci_Cl_cell = 2*(initSetBase.solution.Ci_Cu2_cell+initSetBase.solution.Ci_Fe2_cell)+initSetBase.solution.Ci_H_cell+3*initSetBase.solution.Ci_Fe3_cell;
    initSetBase.solution.Ci_AuCl4_cell = 0.0;
    initSetBase.solution.Ci_cell = [initSetBase.solution.Ci_Cu2_cell initSetBase.solution.Ci_Sn2_cell initSetBase.solution.Ci_Fe2_cell ...
    initSetBase.solution.Ci_Fe3_cell initSetBase.solution.Ci_Ag_cell initSetBase.solution.Ci_Au3_cell ...
    initSetBase.solution.Ci_Pd2_cell initSetBase.solution.Ci_H_cell initSetBase.solution.Ci_Cl_cell initSetBase.solution.Ci_AuCl4_cell];

    %leching vessel concentrations (extraction)
    initSetBase.solution.Ci_Cu2_lch = initSetBase.solution.Ci_Cu2_cell;
    initSetBase.solution.Ci_Sn2_lch = initSetBase.solution.Ci_Sn2_cell;
    initSetBase.solution.Ci_Fe2_lch = initSetBase.solution.Ci_Fe2_cell;
    initSetBase.solution.Ci_Fe3_lch = initSetBase.solution.Ci_Fe3_cell;
    initSetBase.solution.Ci_Ag_lch = initSetBase.solution.Ci_Ag_cell;
    initSetBase.solution.Ci_Au3_lch = initSetBase.solution.Ci_Au3_cell;
    initSetBase.solution.Ci_Pd2_lch = initSetBase.solution.Ci_Pd2_cell;
    initSetBase.solution.Ci_H_lch = initSetBase.solution.Ci_H_cell;
    initSetBase.solution.Ci_Cl_lch = 2*(initSetBase.solution.Ci_Cu2_lch+initSetBase.solution.Ci_Fe2_lch)+initSetBase.solution.Ci_H_lch+3*initSetBase.solution.Ci_Fe3_lch;
    initSetBase.solution.Ci_AuCl4_lch = initSetBase.solution.Ci_AuCl4_cell;
    initSetBase.solution.Ci_lch = [initSetBase.solution.Ci_Cu2_lch initSetBase.solution.Ci_Sn2_lch initSetBase.solution.Ci_Fe2_lch ... 
    initSetBase.solution.Ci_Fe3_lch initSetBase.solution.Ci_Ag_lch initSetBase.solution.Ci_Au3_lch ...
    initSetBase.solution.Ci_Pd2_lch initSetBase.solution.Ci_H_lch initSetBase.solution.Ci_Cl_lch initSetBase.solution.Ci_AuCl4_lch];
    %Base metal setup
    paramSetBase = struct;
    paramSetBase.temp = 298; %K
    paramSetBase.pres = 1; % atm
    paramSetBase.Q = randintF(2,10,0);% L/s (flowrate)
    %cell dimension information
    paramSetBase.length = randintF(3,10,0); % m length of electrodes in flow direction x
    paramSetBase.height = randintF(0.5,2,0); % m height of electrodes    paramSetBase.spacing_x = 0.1; % m gap between end of electrode and vessel inlet/outlet
    paramSetBase.spacing_y = 0.045; %m spacing between electrodes 
    paramSetBase.spacing_x = 0.05; %m spacing between end of cell and electrodes
    paramSetBase.n_units = randintF(2,8,1); %number of anode-cathode surface pairs
    paramSetBase.vol_cell = (paramSetBase.n_units*paramSetBase.spacing_y*...
    paramSetBase.height*(paramSetBase.length+2*paramSetBase.spacing_x))*1000; %L Volume of electrolyte in cell
    %Electrode areas (one side), cm^2
    paramSetBase.S_cat = (paramSetBase.height*100)*(paramSetBase.length*100);
    paramSetBase.S_an = paramSetBase.S_cat;
    %Cross sectional area of cell
    paramSetBase.A_cell = paramSetBase.S_cat;
    %L (Initial) volume of bed holding the particles assuming the bed is half
    %full
    vol_bed = (V_PCB_total/0.6/0.5);
    paramSetBase.vol_lch = vol_bed-V_PCB_total; %L, volume of electrolyte in bed

    paramSetBase.mode = 1; %1 - potentiostat, 2 - galvanostat
    %Applied Voltage (potentiostat)
    paramSetBase.V_app = randintF(2,6,0); %V
    %Applied Current to Cell (Galvanostat)
    paramSetBase.I_app = 36*0.01414; %A
    %Processing time
    paramSetBase.tfinal = tfinal; %s

    %Max current density for all rxns
    paramSetBase.iL_default = 1; %A*m/dm^3
    %fsolve options
    paramSetBase.foptions = optimoptions(@fsolve, 'Display','off', ...
    'MaxFunctionEvaluations', 5000, 'Algorithm', 'trust-region-dogleg', 'StepTolerance', 1E-7);

    %cathode calculations
    thicc_cat = 0.01; %m, thickness of cathode
    A_cat = paramSetBase.S_cat/2*10^-4;  %area of one side in m^2
    s = A_cat^0.5; %assuming square shape
    V_cat = A_cat*thicc_cat; %cathode volume in m^3
    m_cat = rho(6)*V_cat; %Iron mass
    initSetBase.m_deposited = [0 0 m_cat 0 0 0]; %Cu Sn Fe

    try
        resultsBase = metalER(initSetBase,paramSetBase);
        if resultsBase.t(end)<(tfinal-100)
            ie = resultsBase.ie(end);
            if ie == 14
                error('Simulation exited early. Complex number detected in last timestep')
            elseif ie == 15
                error('Simulation exited early. ODE Solver took too long.')
            else
                error('Simulation exited early. Unknown reason.')
            end
        
        end
        disp('Base Sim Success');
        SuccessValsBase(run) = 1;
        save(strcat("SuccessBase\Run",datestr(clock,'mmddHHMMSS'),".mat"),'resultsBase')
    catch exception
        disp('Run Failed');
        SuccessValsBase(run) = 0;
        save(strcat("FailBase\Run",datestr(clock,'mmddHHMMSS'),".mat"),'paramSetBase','initSetBase','exception')
    end
    
    %% Precious metal

    solution = 2; %1 is Cl- base metal, 2 is S2O3 precious metal
    propertiesMetals;
    initSetPrecious = struct;
    %characteristics of solid PCB input
    initSetPrecious.solidPCB.r_particles = resultsBase.PCB.r_particles(size(resultsBase.t,1));
    m_PCB_p = resultsBase.PCB.massRem(size(resultsBase.t,1),:);
    initSetPrecious.solidPCB.m_PCB_total = sum(m_PCB_p); %assumed waste input of 100000 kg/yr, 
    initSetPrecious.solidPCB.wtfrac_PCB = m_PCB_p/sum(m_PCB_p);
    initSetPrecious.m_deposited = [eps 0 0 0 0 0]; %Given a small mass of iron to prevent error
    V_PCB_total = sum(initSetPrecious.solidPCB.m_PCB_total.*initSetPrecious.solidPCB.wtfrac_PCB./rho)*1000;%L

    %characteristics of starting solution
    initSetPrecious.solution.type = solution;%1 is Cl- base metal, 2 is S2O3 precious metal
    %initial concentrations in mol/L
    %Cell Concentrations (recovery)
    initSetPrecious.solution.Ci_Cu2_cell = 0.0;
    initSetPrecious.solution.Ci_Sn2_cell = 0.0;
    initSetPrecious.solution.Ci_Fe2_cell = 0.1;
    initSetPrecious.solution.Ci_Fe3_cell = randintF(0.4,1,0);
    initSetPrecious.solution.Ci_Ag_cell = 0.0;
    initSetPrecious.solution.Ci_Au3_cell = 0.0;
    initSetPrecious.solution.Ci_Pd2_cell = 0.0;
    initSetPrecious.solution.Ci_H_cell = 1E-10;
    initSetPrecious.solution.Ci_S2O3_cell = 0.5;
    initSetPrecious.solution.Ci_AuCl4_cell = 0;
    initSetPrecious.solution.Ci_cell = [initSetPrecious.solution.Ci_Cu2_cell initSetPrecious.solution.Ci_Sn2_cell initSetPrecious.solution.Ci_Fe2_cell ...
    initSetPrecious.solution.Ci_Fe3_cell initSetPrecious.solution.Ci_Ag_cell initSetPrecious.solution.Ci_Au3_cell initSetPrecious.solution.Ci_Pd2_cell ...
    initSetPrecious.solution.Ci_H_cell initSetPrecious.solution.Ci_S2O3_cell initSetPrecious.solution.Ci_AuCl4_cell];

    %leching vessel concentrations (extraction)
    initSetPrecious.solution.Ci_Cu2_lch = initSetPrecious.solution.Ci_Cu2_cell;
    initSetPrecious.solution.Ci_Sn2_lch = initSetPrecious.solution.Ci_Sn2_cell;
    initSetPrecious.solution.Ci_Fe2_lch = initSetPrecious.solution.Ci_Fe2_cell;
    initSetPrecious.solution.Ci_Fe3_lch = initSetPrecious.solution.Ci_Fe3_cell;
    initSetPrecious.solution.Ci_Ag_lch = initSetPrecious.solution.Ci_Ag_cell;
    initSetPrecious.solution.Ci_Au3_lch = initSetPrecious.solution.Ci_Au3_cell;
    initSetPrecious.solution.Ci_Pd2_lch = initSetPrecious.solution.Ci_Pd2_cell;
    initSetPrecious.solution.Ci_H_lch = initSetPrecious.solution.Ci_H_cell;
    initSetPrecious.solution.Ci_S2O3_lch = initSetPrecious.solution.Ci_S2O3_cell;
    initSetPrecious.solution.Ci_AuCl4_lch = initSetPrecious.solution.Ci_AuCl4_cell;
    initSetPrecious.solution.Ci_lch = [initSetPrecious.solution.Ci_Cu2_lch initSetPrecious.solution.Ci_Sn2_lch initSetPrecious.solution.Ci_Fe2_lch ... 
    initSetPrecious.solution.Ci_Fe3_lch initSetPrecious.solution.Ci_Ag_lch initSetPrecious.solution.Ci_Au3_lch initSetPrecious.solution.Ci_Pd2_lch ...
    initSetPrecious.solution.Ci_H_lch initSetPrecious.solution.Ci_S2O3_lch initSetPrecious.solution.Ci_AuCl4_lch];
    paramSetPrecious = struct;
    paramSetPrecious.temp = 298; %K
    paramSetPrecious.pres = 1; % atm
    paramSetPrecious.Q = randintF(2,10,0);% L/s (flowrate)
    %cell dimension information
    paramSetPrecious.length = randintF(3,10,0); % m length of electrodes in flow direction x
    paramSetPrecious.height = randintF(0.5,2,0); % m height of electrodes
    paramSetPrecious.spacing_x = 0.1; % m gap between end of electrode and vessel inlet/outlet
    paramSetPrecious.spacing_y = 0.045; %m spacing between electrodes 
    paramSetPrecious.n_units = randintF(2,8,1); %number of anode-cathode surface pairs
    paramSetPrecious.vol_cell = (paramSetPrecious.n_units*paramSetPrecious.spacing_y*...
    paramSetPrecious.height*(paramSetPrecious.length+2*paramSetPrecious.spacing_x))*1000; %L
    %Electrode areas, cm^2
    paramSetPrecious.S_cat = (paramSetPrecious.height*100)*(paramSetPrecious.length*100);
    paramSetPrecious.S_an = paramSetPrecious.S_cat;
    %Cross sectional area of cell
    paramSetPrecious.A_cell = paramSetPrecious.S_cat;
    %L (Initial) volume of bed holding the particles assuming the bed is half
    %full
    vol_bed = (V_PCB_total/0.6/0.5);
    paramSetPrecious.vol_lch = vol_bed-V_PCB_total; %L, volume of electrolyte in bed


    paramSetPrecious.mode = 1; %1 - potentiostat, 2 - galvanostat
    %Applied Voltage (potentiostat)
    paramSetPrecious.V_app = randintF(3,9,0); %V
    %Applied Current to Cell (Galvanostat)
    paramSetPrecious.I_app = 25;%36*0.01414; %A
    paramSetPrecious.tfinal = tfinal; %s

    %Max current density for all rxns
    paramSetPrecious.iL_default = 1; %A/cm^2
    %fsolve options
    paramSetPrecious.foptions = optimoptions(@fsolve, 'Display','off', ...
    'MaxFunctionEvaluations', 5000, 'Algorithm', 'trust-region-dogleg', 'StepTolerance', 1E-7);

    ie = 0;
    try
        resultsPrecious = metalER(initSetPrecious,paramSetPrecious);
        if resultsPrecious.t(end)<(tfinal-100)
            ie = resultsPrecious.ie(numel(resultsPrecious.ie));
            if ie == 14
                error('Simulation exited early. Complex number detected in last timestep')
            elseif ie == 15
                error('Simulation exited early. ODE Solver took too long.')
            else
                error('Simulation exited early. Unknown reason.')
            end
        
        end
        disp('Sim Success');
        SuccessValsPrecious(run) = 1;
        save(strcat("SuccessPrecious\Run",datestr(clock,'mmddHHMMSS'),".mat"),'resultsPrecious')
    catch exception
        disp('Run Failed');
        SuccessValsPrecious(run) = 0;
        save(strcat("FailPrecious\Run",datestr(clock,'mmddHHMMSS'),".mat"),'paramSetPrecious','initSetPrecious','exception')
    end
end
toc

function randNo = randintF(a,b,isInteger)
    % b must always be greater than a
    if isInteger
        randNo = randi(b-a)+a;
    else
        randNo = (b-a)*rand()+a;
    end
end
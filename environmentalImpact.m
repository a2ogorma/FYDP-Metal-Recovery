function [resultsEnvironmental] = environmentalImpact(resultsPreprocessing, resultsBase, resultsPrecious)
    %energyIntensity
    PF = 0.95; %power factor for energy consumed to produced in electrical grid
    resultsEnvironmental.energy.pumps = 8760*resultsPreprocessing.CF*(resultsBase.numberUnits*resultsBase.practical.pump.BHP + resultsPrecious.numberUnits*resultsPrecious.practical.pump.BHP); %make sure in kWh
    resultsEnvironmental.energy.ESP = 8760*resultsPreprocessing.CF*resultsPreprocessing.workingFactor*resultsPreprocessing.ESP.power; %make sure in kWh
    resultsEnvironmental.energy.grinder = 8760*resultsPreprocessing.CF*resultsPreprocessing.workingFactor*resultsPreprocessing.grinder.power; %make sure in kWh
    resultsEnvironmental.energy.drum = 8760*resultsPreprocessing.CF*resultsPreprocessing.workingFactor*resultsPreprocessing.drum.power; %make sure in kWh
    resultsEnvironmental.energy.EWbase = trapz(resultsBase.t,(resultsBase.electrowinning.V_calc*sum(resultsBase.electrowinning.I_an,2)))/(PF*1000*3600); %in kWh
    resultsEnvironmental.energy.EWprecious = trapz(resultsPrecious.t,(resultsPrecious.electrowinning.V_calc*sum(resultsPrecious.electrowinning.I_an,2)))/(PF*1000*3600); %in kWh
    resultsEnvironmental.energy.total = (resultsEnvironmental.energy.pumps + resultsEnvironmental.energy.ESP + resultsEnvironmental.energy.grinder + resultsEnvironmental.energy.drum + resultsEnvironmental.energy.EWbase + resultsEnvironmental.energy.EWprecious);
    resultsEnvironmental.metrics.energyIntensity = resultsEnvironmental.energy.total/(resultsPreprocessing.Throughput/1000); %in kWh/tonne
    %carbonIntensity
    carbonIntensityGrid = 30; %g CO2 eq/kWh from NIR (mentioned in report)
    resultsEnvironmental.metrics.carbonIntensity = resultsEnvironmental.metrics.energyIntensity*carbonIntensityGrid; %in gCO2e/tonne
    %waterIntensity
    cycleWaterBase = 3; %amount of cycles the solution will leach metals through until the solution is drained and replaced by a fresh solution
    cycleWaterPrecious = 3; %amount of cycles the solution will leach metals through until the solution is drained and replaced by a fresh solution
    resultsEnvironmental.water.BaseVolume = 1.05*(resultsBase.init.paramSet.vol_cell + resultsBase.init.paramSet.vol_lch)*resultsBase.numberUnits; %volume at any given time
    resultsEnvironmental.water.BaseCycles = resultsPreprocessing.CF*8760*3600/resultsBase.init.paramSet.tfinal; %total cycles the system go through per year
    resultsEnvironmental.water.PreciousVolume = 1.05*(resultsPrecious.init.paramSet.vol_cell + resultsPrecious.init.paramSet.vol_lch)*resultsPrecious.numberUnits; %volume at any given time
    resultsEnvironmental.water.PreciousCycles = resultsPreprocessing.CF*8760*3600/resultsPrecious.init.paramSet.tfinal; %total cycles the system go through per year
    resultsEnvironmental.water.annualWaterBase = resultsEnvironmental.water.BaseVolume*(resultsEnvironmental.water.BaseCycles/cycleWaterBase);
    resultsEnvironmental.water.annualWaterPrecious = resultsEnvironmental.water.PreciousVolume*(resultsEnvironmental.water.PreciousCycles/cycleWaterPrecious);
    resultsEnvironmental.metrics.waterIntensity = (resultsEnvironmental.water.annualWaterBase+resultsEnvironmental.water.annualWaterPrecious)/(resultsPreprocessing.Throughput/1000); %L/tonne
    %wasteRecovery
    resultsEnvironmental.waste.platedBase = resultsEnvironmental.water.BaseCycles*resultsBase.numberUnits*(sum(resultsBase.electrowinning.m_plated(end,:))-sum(resultsBase.electrowinning.m_plated(1,:)));
    resultsEnvironmental.waste.platedPrecious = resultsEnvironmental.water.PreciousCycles*resultsPrecious.numberUnits*(sum(resultsPrecious.electrowinning.m_plated(end,:))-sum(resultsPrecious.electrowinning.m_plated(1,:)));
    resultsEnvironmental.waste.loadedBase = resultsEnvironmental.water.BaseCycles*resultsBase.numberUnits*resultsBase.init.initSet.solidPCB.m_PCB_total;
    resultsEnvironmental.metrics.wasteRecovery = (resultsEnvironmental.waste.platedBase+resultsEnvironmental.waste.platedPrecious)/resultsEnvironmental.waste.loadedBase;

    %% Economic %%
    CEPCI = 607.5/400;
    CADUSDconv = 1.27;
    %grinder
    x = resultsPreprocessing.grinder.power; %kW
    resultsEconomic.grinder.Cp = 0.0028*x^3-1.254*x^2+269.26*x+10731;
    resultsEconomic.grinder.Fbm = 2.8;
    resultsEconomic.grinder.capcost = CEPCI*CADUSDconv*resultsEconomic.grinder.Cp*resultsEconomic.grinder.Fbm;
    %magneticdrum
    resultsEconomic.drum.Cp = 13000;
    resultsEconomic.drum.Fbm = 1;
    resultsEconomic.drum.capcost = CADUSDconv*resultsEconomic.drum.Cp*resultsEconomic.drum.Fbm;
    %ESP
    resultsEconomic.ESP.Cp = 30000;%random number
    resultsEconomic.ESP.Fbm = 1;
    resultsEconomic.ESP.capcost = CADUSDconv*resultsEconomic.ESP.Cp*resultsEconomic.ESP.Fbm;
    
    % Base Metal stage
    %pumpBase
    resultsEconomic.baseStage.pump = pumpCost(resultsBase.practical.pump.shaftPower,FM,resultsBase.numberUnits,CEPCI,CADUSDconv);
    resultsEconomic.baseStage.drive.Cp = 195.06*resultsBase.practical.pump.BHP^0.9103;
    resultsEconomic.baseStage.drive.Fbm = 1.5;
    resultsEconomic.baseStage.drive.capcost = resultsBase.numberUnits*CEPCI*CADUSDconv*resultsEconomic.baseStage.drive.Cp*resultsEconomic.baseStage.drive.Fbm;
    %vesselBase
    x = resultsBase.init.paramSet.vol_lch;
    A = 188.85;
    B = 0.6677;
    resultsEconomic.baseStage.leaching.Cp = A*x^B;
    resultsEconomic.baseStage.leaching.Fbm = 1.235;
    resultsEconomic.baseStage.leaching.capcost = resultsBase.numberUnits*CEPCI*CADUSDconv*resultsEconomic.baseStage.leaching.Cp*resultsEconomic.baseStage.leaching.Fbm;
    x = resultsBase.init.paramSet.vol_cell;
    resultsEconomic.baseStage.electrowinning.Cp = A*x^B;
    resultsEconomic.baseStage.electrowinning.Fbm = 1.235;
    resultsEconomic.baseStage.electrowinning.capcost = resultsBase.numberUnits*CEPCI*CADUSDconv*resultsEconomic.baseStage.leaching.Cp*resultsEconomic.baseStage.leaching.Fbm;
    %electrodes
    A_cell = resultsBase.init.paramSet.A_cell; %m2, verify this works
    thicknessCathode = 0.05; %m, temporary
    matDensityCathode = 11343;%kg/m3, lead
    matCostCathode = 1.61;%$/kg, lead
    resultsEconomic.baseStage.cathodeCost = resultsBase.numberUnits*(A_cell/2)*thicknessCathode*matDensityCathode*matCostCathode;
    thicknessAnode = 0.05; %m, temporary
    matDensityAnode = 7750; %kg/m3, SS
    matCostAnode = 2.66; %$/kg, SS  
    resultsEconomic.baseStage.anodeCost = resultsBase.numberUnits*(((resultsBase.init.paramSet.n_epairs/2)+1)/(resultsBase.init.paramSet.n_epairs/2))*(A_cell/2)*thicknessAnode*matDensityAnode*matCostAnode;
   
    % Precious Metal stage
    %pumpPrecious
    resultsEconomic.preciousStage.pump = pumpCost(resultsPrecious.practical.pump.shaftPower,FM,resultsPrecious.numberUnits,CEPCI,CADUSDconv);
    resultsEconomic.preciousStage.drive.Cp = 195.06*resultsPrecious.practical.pump.BHP^0.9103;
    resultsEconomic.preciousStage.drive.Fbm = 1.5;
    resultsEconomic.preciousStage.drive.capcost = resultsPrecious.numberUnits*CEPCI*CADUSDconv*resultsEconomic.preciousStage.drive.Cp*resultsEconomic.preciousStage.drive.Fbm;
    %vesselPrecious
    x = resultsPrecious.init.paramSet.vol_lch;
    A = 188.85;
    B = 0.6677;
    resultsEconomic.preciousStage.leaching.Cp = A*x^B;
    resultsEconomic.preciousStage.leaching.Fbm = 1.235;
    resultsEconomic.preciousStage.leaching.capcost = CEPCI*CADUSDconv*resultsEconomic.preciousStage.leaching.Cp*resultsEconomic.preciousStage.leaching.Fbm;
    x = resultsPrecious.init.paramSet.vol_cell;
    resultsEconomic.preciousStage.electrowinning.Cp = A*x^B;
    resultsEconomic.preciousStage.electrowinning.Fbm = 1.235;
    resultsEconomic.preciousStage.electrowinning.capcost = resultsPrecious.numberUnits*CEPCI*CADUSDconv*resultsEconomic.preciousStage.leaching.Cp*resultsEconomic.preciousStage.leaching.Fbm;
    %electrodes
    A_cell = resultsPrecious.init.paramSet.A_cell; %m2, verify this works
    thicknessCathode = 0.05; %m, temporary
    matDensityCathode = 11343;%kg/m3, lead
    matCostCathode = 1.61;%$/kg, lead
    resultsEconomic.preciousStage.cathodeCost = resultsPrecious.numberUnits*(A_cell/2)*thicknessCathode*matDensityCathode*matCostCathode;
    thicknessAnode = 0.05; %m, temporary
    matDensityAnode = 7750; %kg/m3, SS
    matCostAnode = 2.66; %$/kg, SS  
    resultsEconomic.preciousStage.anodeCost = resultsPrecious.numberUnits*(((resultsPrecious.init.paramSet.n_epairs/2)+1)/(resultsPrecious.init.paramSet.n_epairs/2))*(A_cell/2)*thicknessAnode*matDensityAnode*matCostAnode;
    
    %%%%%ADD rectifier costs
    %%%%% add total capital cost amount summation
    
    % Operating
    %Chemical costs
    %solution base
    solutionCost = 09182309183908; %$/L obv change im sleepy. Note that this
    ... should be all components in solution, including any additives and stuff
    resultsEconomic.operating.solutionBaseCost = resultsEnvironmental.water.annualWaterBase*solutionCost;
    %solution precious
    solutionCost = 09182309183908; %$/L obv change im sleepy. Note that this
    ... should be all components in solution, including any additives and stuff
    resultsEconomic.operating.solutionPreciousCost = resultsEnvironmental.water.annualWaterPrecious*solutionCost;
    %sodium hydroxide addition for pH balancing?
    
    %Labour costs
    labourCost = 17; %$/hr
    %add more labour stuff here
    
    %Electrical
    ElecCost = 0.139; % $/kWh
    resultsEconomic.operating.energyCost =  resultsEnvironmental.energy.total*ElecCost;
    
    %overall operating
    
    % Moneymaking
    %deposited metal selling amounts
    
    % Economics finishing metrics
    
end

function [pump] = pumpCost(power,FM,numberUnits,CEPCI,CADUSDconv)
    %power in kW
    pump.Cp = 4240*power^0.3624;
    pump.FM = FM;
    pump.Fbm = 1.4957*FM+2.1287;
    pump.capcost = numberUnits*CEPCI*CADUSDconv*pump.Cp*pump.Fbm;
end
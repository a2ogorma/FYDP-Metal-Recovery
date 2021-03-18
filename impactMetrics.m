function [resultsEnvironmental, resultsEconomic] = impactMetrics(resultsPreprocessing, resultsBase, resultsPrecious)
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
    resultsEconomic.preprocessing.grinder.Cp = 0.0028*x^3-1.254*x^2+269.26*x+10731;
    resultsEconomic.preprocessing.grinder.Fbm = 2.8;
    resultsEconomic.preprocessing.grinder.capcost = CEPCI*CADUSDconv*resultsEconomic.preprocessing.grinder.Cp*resultsEconomic.preprocessing.grinder.Fbm;
    %magneticdrum
    resultsEconomic.preprocessing.drum.Cp = 13000;
    resultsEconomic.preprocessing.drum.Fbm = 1;
    resultsEconomic.preprocessing.drum.capcost = CADUSDconv*resultsEconomic.preprocessing.drum.Cp*resultsEconomic.preprocessing.drum.Fbm;
    %ESP
    resultsEconomic.preprocessing.ESP.Cp = 30000;%random number
    resultsEconomic.preprocessing.ESP.Fbm = 1;
    resultsEconomic.preprocessing.ESP.capcost = CADUSDconv*resultsEconomic.preprocessing.ESP.Cp*resultsEconomic.preprocessing.ESP.Fbm;
    %preprocessingtotal
    resultsEconomic.preprocessing.capcost = resultsEconomic.preprocessing.grinder.capcost + resultsEconomic.preprocessing.drum.capcost + resultsEconomic.preprocessing.ESP.capcost;
    
    % Base Metal stage
    %pumpBase
    FM = 1;
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
    resultsEconomic.baseStage.anodeCost = resultsBase.numberUnits*(((resultsBase.init.paramSet.n_units/2)+1)/(resultsBase.init.paramSet.n_units/2))*(A_cell/2)*thicknessAnode*matDensityAnode*matCostAnode;
    %basecapcosttotal
    resultsEconomic.baseStage.capcost = resultsEconomic.baseStage.pump.capcost + resultsEconomic.baseStage.drive.capcost + resultsEconomic.baseStage.leaching.capcost + resultsEconomic.baseStage.electrowinning.capcost + resultsEconomic.baseStage.cathodeCost + resultsEconomic.baseStage.anodeCost;
    
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
    resultsEconomic.preciousStage.cathodecost = resultsPrecious.numberUnits*(A_cell/2)*thicknessCathode*matDensityCathode*matCostCathode;
    thicknessAnode = 0.05; %m, temporary
    matDensityAnode = 7750; %kg/m3, SS
    matCostAnode = 2.66; %$/kg, SS  
    resultsEconomic.preciousStage.anodecost = resultsPrecious.numberUnits*(((resultsPrecious.init.paramSet.n_units/2)+1)/(resultsPrecious.init.paramSet.n_units/2))*(A_cell/2)*thicknessAnode*matDensityAnode*matCostAnode;
    %preciouscapcosttotal
    resultsEconomic.preciousStage.capcost = resultsEconomic.preciousStage.pump.capcost + resultsEconomic.preciousStage.drive.capcost + resultsEconomic.preciousStage.leaching.capcost + resultsEconomic.preciousStage.electrowinning.capcost + resultsEconomic.preciousStage.cathodecost + resultsEconomic.preciousStage.anodecost;
    
    %%%%%ADD rectifier costs maybe?
    %totalcapcost
    resultsEconomic.capcost = resultsEconomic.preciousStage.capcost + resultsEconomic.baseStage.capcost + resultsEconomic.preprocessing.capcost;
    contingencyandfees = 0.18;
    greenfield = 0.3; %auxiliary facilities and utilities supply
    resultsEconomic.fixedCost = (1+contingencyandfees+greenfield)*resultsEconomic.capcost;
    workingCapital = 0.15; % 10-20%
    resultsEconomic.workingCapital = workingCapital*resultsEconomic.fixedCost;
    resultsEconomic.totalCapitalInvestment = resultsEconomic.fixedCost+resultsEconomic.workingCapital;
    
    % Operating
    %Chemical costs
    %solution base
    solutionCost = 5; %$/L obv change im sleepy. Note that this
    ... should be all components in solution, including any additives and stuff
    resultsEconomic.operating.solutionBaseCost = resultsEnvironmental.water.annualWaterBase*solutionCost;
    %solution precious
    solutionCost = 5; %$/L obv change im sleepy. Note that this
    ... should be all components in solution, including any additives and stuff
    resultsEconomic.operating.solutionPreciousCost = resultsEnvironmental.water.annualWaterPrecious*solutionCost;
    %sodium hydroxide addition for pH balancing?
    
    %Labour costs
    labourCost = 20; %$/hr - suggestion from 480
    labourDensity = 3;%people working while plant operating
    labourHours = 8760*resultsPreprocessing.CF*labourDensity; %Assumed 3 people working per working hours rn
    resultsEconomic.operating.labourCost = labourHours*labourCost;
    resultsEconomic.operating.supervisorCost = 0.15*resultsEconomic.operating.labourCost; %10-30% supervisory and clerical according to 480
    labCharges = 0.15;
    resultsEconomic.operating.labCharges = labCharges*resultsEconomic.operating.labourCost;
    %Electrical
    ElecCost = 0.139; % $/kWh
    resultsEconomic.operating.energyCost =  resultsEnvironmental.energy.total*ElecCost;
    
    %Waste disposal
    wastewaterDisposal = 1.3;%dummy val, $/L
    resultsEconomic.operating.wastewaterDisposal = wastewaterDisposal*(resultsEnvironmental.water.annualWaterPrecious+resultsEnvironmental.water.annualWaterBase);
    %Other direct expenses
    mainrepair = 0.04; %2-10% of fixed capital for maintenance and repair
    resultsEconomic.operating.mainrepair = mainrepair*resultsEconomic.fixedCost;
    opersupplies = 0.15; %10-20% of maintenance and repair
    resultsEconomic.operating.opersupplies = opersupplies*resultsEconomic.operating.mainrepair;
    resultsEconomic.operating.subtotal = resultsEconomic.operating.solutionBaseCost + resultsEconomic.operating.solutionPreciousCost + resultsEconomic.operating.labourCost + resultsEconomic.operating.supervisorCost + resultsEconomic.operating.labCharges + resultsEconomic.operating.energyCost + resultsEconomic.operating.wastewaterDisposal + resultsEconomic.operating.mainrepair + resultsEconomic.operating.opersupplies;
    patentsRoyalties = 0.03; % 0-6% of total expense
    resultsEconomic.operating.patentsRoyalties = resultsEconomic.operating.subtotal*patentsRoyalties;
    resultsEconomic.operating.totalDirectExpenses = resultsEconomic.operating.subtotal + resultsEconomic.operating.patentsRoyalties;
    %indirect expenses
    overheadRate = 0.6; %50-70% of oplabour, supervision, maintenance
    resultsEconomic.operating.overhead = overheadRate*(resultsEconomic.operating.labourCost + resultsEconomic.operating.supervisorCost + resultsEconomic.operating.mainrepair);
    localTaxes = 0.02; %1-3 of fixed cost
    resultsEconomic.operating.localTaxes = localTaxes*resultsEconomic.fixedCost;
    insurance = 0.01; % 1-2 of fixed cost
    resultsEconomic.operating.insurance = insurance*resultsEconomic.fixedCost;
    resultsEconomic.totalIndirectExpenses = resultsEconomic.operating.overhead + resultsEconomic.operating.localTaxes + resultsEconomic.operating.insurance;
    %total MFG
    resultsEconomic.totalManufacturingExpenses = resultsEconomic.totalIndirectExpenses +  resultsEconomic.operating.totalDirectExpenses;
    
    admin = 0.25; % of overhead
    resultsEconomic.operating.admin = admin*resultsEconomic.operating.overhead;
    distSelling = 0.10; % of total expense
    resultsEconomic.operating.distSelling = distSelling*resultsEconomic.totalManufacturingExpenses;
    RDcosts = 0.05;% of total expense
    resultsEconomic.operating.RDcosts = RDcosts*resultsEconomic.totalManufacturingExpenses;
    depreciation = 0.10; %10 of fixed capital
    resultsEconomic.operating.depreciation = depreciation*resultsEconomic.fixedCost;
    resultsEconomic.totalGeneralExpenses = resultsEconomic.operating.admin + resultsEconomic.operating.distSelling + resultsEconomic.operating.RDcosts + resultsEconomic.operating.depreciation;
    
    resultsEconomic.totalExpenses = resultsEconomic.totalGeneralExpenses+resultsEconomic.totalManufacturingExpenses;
    
    % Moneymaking
    %deposited metal selling amounts
    metalPrices = [4.77 38.8896 1.7 1039.55 69335.99 95246.27]; %find better sources for this
    sellingRate = 0.9;% 90% of metal price
    platedBase = subplus(resultsBase.electrowinning.m_plated(end,:)-resultsBase.electrowinning.m_plated(1,:));
    resultsEconomic.revenue.Base = resultsEnvironmental.water.BaseCycles*platedBase*metalPrices'*sellingRate;
    platedPrecious = (resultsPrecious.electrowinning.m_plated(end,:)-resultsPrecious.electrowinning.m_plated(1,:));
    resultsEconomic.revenue.Precious = resultsEnvironmental.water.PreciousCycles*platedPrecious*metalPrices'*sellingRate;    
    resultsEconomic.totalRevenue = resultsEconomic.revenue.Base + resultsEconomic.revenue.Precious;
    
    resultsEconomic.netAnnualbeforeTax = resultsEconomic.totalRevenue - resultsEconomic.totalExpenses;
    taxRate = 0.3;
    resultsEconomic.metrics.netAnnualafterTax = (1-taxRate)*resultsEconomic.netAnnualbeforeTax;
    resultsEconomic.metrics.percentOp_of_Rev = resultsEconomic.totalExpenses/resultsEconomic.totalRevenue;
    resultsEconomic.metrics.paybackPeriod = resultsEconomic.fixedCost/resultsEconomic.metrics.netAnnualafterTax;
    %maybe include discount rate for NPV calcs?
end

function [pump] = pumpCost(power,FM,numberUnits,CEPCI,CADUSDconv)
    %power in kW
    pump.Cp = 4240*power^0.3624;
    pump.FM = FM;
    pump.Fbm = 1.4957*FM+2.1287;
    pump.capcost = numberUnits*CEPCI*CADUSDconv*pump.Cp*pump.Fbm;
end
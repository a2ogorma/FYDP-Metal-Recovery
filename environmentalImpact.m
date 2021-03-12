function [resultsEnvironmental] = environmentalImpact(resultsPreprocessing, resultsBase, resultsPrecious)
    %energyIntensity
    PF = 0.95; %power factor for energy consumed to produced in electrical grid
    resultsEnvironmental.energy.pumps = resultsPreprocessing.CF*(resultsBase.numberUnits*resultsBase.practical.pump.shaftPower + resultsPrecious.numberUnits*resultsPrecious.practical.pump.shaftPower); %make sure in kWh
    resultsEnvironmental.energy.ESP = resultsPreprocessing.CF*resultsPreprocessing.workingHours*resultsPreprocessing.ESP.power; %make sure in kWh
    resultsEnvironmental.energy.grinder = resultsPreprocessing.CF*resultsPreprocessing.workingHours*resultsPreprocessing.grinder.power; %make sure in kWh
    resultsEnvironmental.energy.drum = resultsPreprocessing.CF*resultsPreprocessing.workingHours*resultsPreprocessing.drum.power; %make sure in kWh
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
end
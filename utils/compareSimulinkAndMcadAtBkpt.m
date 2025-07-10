function maxTempError = compareSimulinkAndMcadAtBkpt(modelName, mcadIntf, coolingSystemsEnabled, bkptIdxComb, torqueVal)
%   COMPARESIMULINKANDMCADATBKPT
%
%   Works with an arbitrary list of cooling systems supplied in
%   COOLINGSYSTEMSENABLED.   The first element in BKPTIDXCOMB is always the
%   speed index; every cooling branch then contributes *two* indices
%   (flow-rate, inlet-temperature) in the same order as
%   COOLINGSYSTEMSENABLED.
%
%   Copyright 2022-2025 The MathWorks, Inc.

%% --- basic checks -------------------------------------------------------
nSys          = numel(coolingSystemsEnabled);
expectedLen   = 1 + 2*nSys;          % 1 speed  + 2 per system
assert(numel(bkptIdxComb)==expectedLen, ...
    "bkptIdxComb must contain %d indices (1 speed + 2 per cooling system).",expectedLen);

%% --- grab sampling grid from model workspace ---------------------------
mdlWks        = get_param(modelName,'ModelWorkspace');
SSND          = mdlWks.getVariable('StateSpaceND');
SG            = SSND.SamplingGrid;

% helpers
getVec  = @(field) unique(SG.(field)(:))';            % collapse N-D array
mkName  = @(s) regexprep(s,'[^A-Za-z0-9_]','');       % "Housing Water Jacket" → "HousingWaterJacket", "Spray_RadialRotor" → "Spray_RadialRotor", etc.

%% --- map indices → physical values -------------------------------------
rpmIdx         = bkptIdxComb(1);
rpmVec         = getVec('w');
rpmVal         = rpmVec(rpmIdx);

flowVals  = zeros(1,nSys);
tempVals  = zeros(1,nSys);
blkBases  = cell (1,nSys);

for k = 1:nSys
    blkBases{k}   = mkName(coolingSystemsEnabled{k});

    flowVec       = getVec(sprintf('fr%d',k));
    tempVec       = getVec(sprintf('Tin%d',k));

    flowVals(k)   = flowVec(bkptIdxComb(1+2*(k-1)+1));
    tempVals(k)   = tempVec(bkptIdxComb(1+2*(k-1)+2));
end

%% --- feed values into Simulink -----------------------------------------
set_param([modelName '/SpeedRPM'],  'Value',num2str(rpmVal));
set_param([modelName '/TorqueNm'],  'Value',num2str(torqueVal));

for k = 1:nSys
    blk = blkBases{k};
    set_param([modelName '/' blk '_FlowRate_lpm'   ],'Value',num2str(flowVals(k)));
    set_param([modelName '/' blk '_InletTemp_degC' ],'Value',num2str(tempVals(k)));
end

stopTime = 1000;                           % [s]
out      = sim(modelName,'StopTime',num2str(stopTime));

TnodesSeries  = out.yout{1}.Values;        % Simulink temperatures
TnodesInit    = mdlWks.getVariable('TnodesInit');

%% --- configure & run Motor-CAD -----------------------------------------
mcadIntf.Tambient_degC      = TnodesInit(1);
mcadIntf.Shaft_Speed_RPM    = rpmVal;

for k = 1:nSys
    base            = blkBases{k};
    % Most Motor-CAD models use “…FlowRate_m3ps” / “…InletTemperature_degC”.
    % (Falls back to "…FlowVelocity_mps" if needed.)
    propFlow        = sprintf('%s_FlowRate_m3ps',base);
    if isprop(mcadIntf,propFlow)
        mcadIntf.(propFlow) = flowVals(k)/60/1000; % lpm → m³/s
    else
        propFlow = sprintf('%s_FlowVelocity_mps',base);
        mcadIntf.(propFlow) = flowVals(k); % m/s
    end
    propTin = sprintf('%s_InletTemperature_degC',base);
    if isprop(mcadIntf,propTin)
        mcadIntf.(propTin) = tempVals(k);
    else
        % multi-nozzle, with Front and Rear - set both
        propTinF = sprintf('%s_InletTemperature_F_degC',base);
        mcadIntf.(propTinF) = tempVals(k);
        propTinR = sprintf('%s_InletTemperature_R_degC',base);
        mcadIntf.(propTinR) = tempVals(k);
    end
end

mcadIntf.EnableStatorTempCoeffRes = 1;
mcadIntf.EnableRotorTempCoeffRes  = 1;
mcadIntf.updateModel();

numTimeSteps = 50;
mcadIntf.runThermalTransientWithSpecifiedTorqueSpeed( ...
    torqueVal,rpmVal,stopTime,numTimeSteps);

allMcadIdxs   = [mcadIntf.NodeNamesAndMcadIdx{:,2}];
[tVecMcad,Tmcad] = mcadIntf.getTransientTemperatureForNodeMcadIdxs(allMcadIdxs);

%% --- error calculation & visualisation ----------------------------------
TmcadSeries = timeseries(Tmcad',tVecMcad);  TmcadSeries.Name = 'Motor-CAD';
[TsMcad,TsSL] = synchronize(TmcadSeries,TnodesSeries,'Union');

err  = squeeze(TsSL.Data)' - TsMcad.Data;
err(TsSL.Time < stopTime/10,:) = [];        % ignore “warm-up” phase
maxTempError = max(abs(err(:)));

figure;
plot(TsSL ,'b'); hold on;
plot(TsMcad,'r--'); grid on;
legend('Simulink','Motor-CAD','Location','best');
title(composeTitle(rpmVal,flowVals,tempVals,coolingSystemsEnabled));
ylabel('Temperature [°C]'); xlabel('Time [s]');

end  % --------------------------------------------------------------------

%% helper: dynamic plot title
function txt = composeTitle(rpm,flw,tmp,sys)
    parts = ["w = " + rpm + " rpm"];
    for i = 1:numel(sys)
        parts(end+1) = sprintf('%s: %.2g lpm / %.2g °C',sys{i},flw(i),tmp(i));
    end
    txt = strjoin(parts,', ');
end
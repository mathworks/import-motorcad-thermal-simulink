function maxTempError = compareSimulinkAndMcadAtBkpt(modelName, mcadIntf, bkptIdxComb, torqueVal)
    % COMPARESIMULINKANDMCADATBKPT returns the maximum temperature absolute
    % difference between Simulink and Motor-CAD in a transient simulation 
    % with specified torque, at a particular combination of breakpoints.
    % It also plots the two simulations superimposed.
    %
    % Input arguments:
    % - modelName: [string/char]: model name of the SROTM, previously generated
    % with generateSimulinkReducedOrderModel
    % - mcadIntf: [mcadinterface.ThermalInterface object]: pre-loaded ThermalInterface object
    % - bkptIdxComb: [int]: Array of breakpoint indices of speed, flow
    % rate, and inlet temperature to use for the comparison.
    % - torqueVal: [double]: Value of shaft torque for the comparison.
    % Output arguments:
    % - maxTempError: [double]: max(abs(Tsrotm-Tmcad)) where Tsimulink
    % is the SROTM node temperatures and Tmcad is the Motor-CAD node
    % temperature, in the transient simulation.

    % Copyright 2022-2025 The MathWorks, Inc.
    
    mdlWks = get_param(modelName,'ModelWorkspace');
    StateSpaceND = mdlWks.getVariable('StateSpaceND');
    samplingGrid = StateSpaceND.SamplingGrid;

    % Obtain Simulink results ------------------------------------------
    if length(bkptIdxComb)==3 % one cooling system - e5_IM_HWJ case
        wBkpts = squeeze(samplingGrid.w(:,1,1))';
        fr1Bkpts = squeeze(samplingGrid.fr1(1,:,1));
        Tin1Bkpts = squeeze(samplingGrid.Tin1(1,1,:))';
        rpmVal = wBkpts(bkptIdxComb(1));
        fr1Val = fr1Bkpts(bkptIdxComb(2));
        Tin1Val = Tin1Bkpts(bkptIdxComb(3));
        set_param(strcat(modelName, '/HousingWaterJacket_Flowrate_lpm'), 'Value', num2str(fr1Val));
        set_param(strcat(modelName, '/HousingWaterJacket_InletTemp_degC'), 'Value', num2str(Tin1Val));
    elseif length(bkptIdxComb)==5 % two cooling systems - e8_IPMSM_HWJandVent case
        wBkpts = squeeze(samplingGrid.w(:,1,1,1,1))';
        fr1Bkpts = squeeze(samplingGrid.fr1(1,:,1,1,1));
        fr2Bkpts = squeeze(samplingGrid.fr2(1,1,:,1,1))';
        Tin1Bkpts = squeeze(samplingGrid.Tin1(1,1,1,:,1))';
        Tin2Bkpts = squeeze(samplingGrid.Tin2(1,1,1,1,:))';
        rpmVal = wBkpts(bkptIdxComb(1));
        fr1Val = fr1Bkpts(bkptIdxComb(2));       
        fr2Val = fr2Bkpts(bkptIdxComb(3));
        Tin1Val = Tin1Bkpts(bkptIdxComb(4));
        Tin2Val = Tin2Bkpts(bkptIdxComb(5));
        set_param(strcat(modelName, '/HousingWaterJacket_Flowrate_lpm'), 'Value', num2str(fr1Val));
        set_param(strcat(modelName, '/HousingWaterJacket_InletTemp_degC'), 'Value', num2str(Tin1Val));
        set_param(strcat(modelName, '/Ventilated_Flowrate_lpm'), 'Value', num2str(fr2Val));
        set_param(strcat(modelName, '/Ventilated_InletTemp_degC'), 'Value', num2str(Tin2Val));
    elseif length(bkptIdxComb)==7 % three cooling systems - e8_IPMSM_SprayMultiNozzle case
        wBkpts = squeeze(samplingGrid.w(:,1,1,1,1,1,1))';
        fr1Bkpts = squeeze(samplingGrid.fr1(1,:,1,1,1,1,1)); 
        Tin1Bkpts = squeeze(samplingGrid.Tin1(1,1,:,1,1,1,1))'; 
        fr2Bkpts = squeeze(samplingGrid.fr2(1,1,1,:,1,1,1))'; 
        Tin2Bkpts = squeeze(samplingGrid.Tin2(1,1,1,1,:,1,1))'; 
        fr3Bkpts = squeeze(samplingGrid.fr3(1,1,1,1,1,:,1))'; 
        Tin3Bkpts = squeeze(samplingGrid.Tin3(1,1,1,1,1,1,:))'; 
        rpmVal = wBkpts(bkptIdxComb(1));
        fr1Val = fr1Bkpts(bkptIdxComb(2));       
        fr2Val = fr2Bkpts(bkptIdxComb(3));
        fr3Val = fr3Bkpts(bkptIdxComb(4));
        Tin1Val = Tin1Bkpts(bkptIdxComb(5));
        Tin2Val = Tin2Bkpts(bkptIdxComb(6));
        Tin3Val = Tin3Bkpts(bkptIdxComb(7));
        set_param(strcat(modelName, '/Spray_RadialHousing_Flowrate_lpm'), 'Value', num2str(fr1Val));
        set_param(strcat(modelName, '/Spray_RadialHousing_InletTemp_degC'), 'Value', num2str(Tin1Val));
        set_param(strcat(modelName, '/HousingWaterJacket_Flowrate_lpm'), 'Value', num2str(fr2Val));
        set_param(strcat(modelName, '/HousingWaterJacket_InletTemp_degC'), 'Value', num2str(Tin2Val));
        set_param(strcat(modelName, '/Spray_RadialRotor_Flowrate_lpm'), 'Value', num2str(fr3Val));
        set_param(strcat(modelName, '/Spray_RadialRotor_InletTemp_degC'), 'Value', num2str(Tin3Val));
    end

    set_param(strcat(modelName, '/TorqueNm'), 'Value', num2str(torqueVal));
    set_param(strcat(modelName, '/SpeedRPM'), 'Value', num2str(rpmVal));
    
    stopTime = 1000; % s
    out = sim(modelName, 'StopTime', num2str(stopTime));
    TnodesSeries = out.yout{1}.Values;
    TnodesInit = mdlWks.getVariable('TnodesInit');

    % Obtain Motor-CAD baseline ------------------------------------------
    mcadIntf.Tambient_degC = TnodesInit(1);
    mcadIntf.Shaft_Speed_RPM = rpmVal;
    if length(bkptIdxComb)==3 % one cooling system - e5_IM_HWJ case
        mcadIntf.HousingWaterJacket_FlowRate_m3ps = fr1Val/60/1000; % lpm to m3ps
        mcadIntf.HousingWaterJacket_InletTemperature_degC = Tin1Val;
    elseif length(bkptIdxComb)==5 % two cooling systems - e8_IPMSM_HWJandVent case
        mcadIntf.HousingWaterJacket_FlowRate_m3ps = fr1Val/60/1000; % lpm to m3ps
        mcadIntf.HousingWaterJacket_InletTemperature_degC = Tin1Val;
        mcadIntf.Ventilated_FlowRate_m3ps = fr2Val/60/1000; % lpm to m3ps
        mcadIntf.Ventilated_InletTemperature_degC = Tin2Val;
    elseif length(bkptIdxComb)==7 % three cooling systems - e8_IPMSM_SprayMultiNozzle case
        mcadIntf.Spray_RadialHousing_FlowRate_m3ps = fr1Val/60/1000; % lpm to m3ps
        mcadIntf.Spray_RadialHousing_InletTemperature_F_degC = Tin1Val;
        mcadIntf.Spray_RadialHousing_InletTemperature_R_degC = Tin1Val;
        mcadIntf.HousingWaterJacket_FlowRate_m3ps = fr2Val/60/1000; % lpm to m3ps
        mcadIntf.HousingWaterJacket_InletTemperature_degC = Tin2Val;
        mcadIntf.Spray_RadialRotor_FlowRate_m3ps = fr3Val/60/1000; % lpm to m3ps
        mcadIntf.Spray_RadialRotor_InletTemperature_F_degC = Tin3Val;
        mcadIntf.Spray_RadialRotor_InletTemperature_R_degC = Tin3Val;
    end

    mcadIntf.EnableStatorTempCoeffRes = 1;
    mcadIntf.EnableRotorTempCoeffRes = 1;
    mcadIntf.updateModel();    
    numTimeSteps = 50;
    mcadIntf.runThermalTransientWithSpecifiedTorqueSpeed(torqueVal, rpmVal, stopTime, numTimeSteps)
    allMcadIdxs = [mcadIntf.NodeNamesAndMcadIdx{:,2}];
    [tVecMcad, TnodesMcad] = mcadIntf.getTransientTemperatureForNodeMcadIdxs(allMcadIdxs);

    % Plot results ------------------------------------------------------
    TnodesMcadSeries = timeseries(TnodesMcad', tVecMcad);
    [TnodesMcadSeries,TnodesSeries] = synchronize(TnodesMcadSeries,TnodesSeries,'Union');
    TnodesMcadSeries.Name = 'Motor-CAD';
    TnodesSeries.Name = 'Simulink';

    TnodesError = squeeze(TnodesSeries.Data)'-TnodesMcadSeries.Data;  
    TnodesError(TnodesSeries.Time < stopTime/10,:) = []; % Ignore initial transient i.e. initial 1/10th of simulation
    maxTempError = max(abs(TnodesError(:)));

    figure();
    h1 = plot(TnodesSeries, 'b');
    hold on
    h2 = plot(TnodesMcadSeries, 'r--');
    hold off
    legend([h1(1), h2(1)], {'Simulink', 'Motor-CAD'});
    if length(bkptIdxComb)==3
        title(strcat('w = ', num2str(rpmVal), ...
            ' rpm, fr = ', num2str(fr1Val), ' lpm, Tin = ', num2str(Tin1Val), ' degC'));
    elseif length(bkptIdxComb)==5
        title(strcat('w = ', num2str(rpmVal), ...
            ' rpm, fr1 = ', num2str(fr1Val), ' lpm, fr2 = ', num2str(fr2Val), ...
            ' lpm, Tin1 = ', num2str(Tin1Val), ' degC, Tin2 = ', num2str(Tin2Val), ' degC'));
    end
    ylabel('Node temperatures [degC]');
    xlabel('Time [s]'); 
    grid on
    
end


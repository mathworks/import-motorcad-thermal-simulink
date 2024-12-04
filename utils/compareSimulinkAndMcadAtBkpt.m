function maxTempError = compareSimulinkAndMcadAtBkpt(modelName, motFile, bkptIdxComb, torqueVal)
    % COMPARESIMULINKANDMCADATBKPT returns the maximum temperature absolute
    % difference between Simulink and Motor-CAD in a transient simulation 
    % with specified torque, at a particular combination of breakpoints.
    % It also plots the two simulations superimposed.
    %
    % Input arguments:
    % - modelName: [string/char]: model name of the SROTM, previously generated
    % with generateSimulinkReducedOrderModel
    % - motFile: [string/char]: Motor-CAD .mot file to compare against.
    % - bkptIdxComb: [int]: Array of breakpoint indices of speed, flow
    % rate, and inlet temperature to use for the comparison.
    % - torqueVal: [double]: Value of shaft torque for the comparison.
    % Output arguments:
    % - maxTempError: [double]: max(abs(Tsrotm-Tmcad)) where Tsimulink
    % is the SROTM node temperatures and Tmcad is the Motor-CAD node
    % temperature, in the transient simulation.

    % Copyright 2022 The MathWorks, Inc.
    
    mdlWks = get_param(modelName,'ModelWorkspace');
    StateSpaceND = mdlWks.getVariable('StateSpaceND');
    samplingGrid = StateSpaceND.SamplingGrid;

    if length(bkptIdxComb)==3 % one cooling system
        wBkpts = squeeze(samplingGrid.w(:,1,1))';
        fr1Bkpts = squeeze(samplingGrid.fr1(1,:,1));
        Tin1Bkpts = squeeze(samplingGrid.Tin1(1,1,:))';
        rpmVal = wBkpts(bkptIdxComb(1));
        fr1Val = fr1Bkpts(bkptIdxComb(2));
        Tin1Val = Tin1Bkpts(bkptIdxComb(3));        
    elseif length(bkptIdxComb)==5 % two cooling systems
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
    end

    set_param(strcat(modelName, '/TorqueNm'), 'Value', num2str(torqueVal));
    set_param(strcat(modelName, '/SpeedRPM'), 'Value', num2str(rpmVal));
    set_param(strcat(modelName, '/HousingWaterJacket_Flowrate_lpm'), 'Value', num2str(fr1Val));
    set_param(strcat(modelName, '/HousingWaterJacket_InletTemp_degC'), 'Value', num2str(Tin1Val));
    if length(bkptIdxComb)==5
        set_param(strcat(modelName, '/Ventilated_Flowrate_lpm'), 'Value', num2str(fr2Val));
        set_param(strcat(modelName, '/Ventilated_InletTemp_degC'), 'Value', num2str(Tin2Val));
    end
    
    stopTime = 1000; % s
    out = sim(modelName, 'StopTime', num2str(stopTime));
    TnodesSeries = out.yout{1}.Values;

    McadIntf = mcadinterface.ThermalInterface(motFile);
    TnodesInit = mdlWks.getVariable('TnodesInit');
    McadIntf.Tambient_degC = TnodesInit(1);
    McadIntf.Shaft_Speed_RPM = rpmVal;
    McadIntf.HousingWaterJacket_FlowRate_m3ps = fr1Val/60/1000; % lpm to m3ps
    McadIntf.HousingWaterJacket_InletTemperature_degC = Tin1Val;
    if length(bkptIdxComb)==5
        McadIntf.Ventilated_FlowRate_m3ps = fr2Val/60/1000; % lpm to m3ps
        McadIntf.Ventilated_InletTemperature_degC = Tin2Val;
    end
    McadIntf.EnableStatorTempCoeffRes = 1;
    McadIntf.EnableRotorTempCoeffRes = 1;
    McadIntf.updateModel();    
    numTimeSteps = 50;
    McadIntf.runThermalTransientWithSpecifiedTorqueSpeed(torqueVal, rpmVal, stopTime, numTimeSteps)
    allMcadIdxs = [McadIntf.NodeNamesAndMcadIdx{:,2}];
    [tVecMcad, TnodesMcad] = McadIntf.getTransientTemperatureForNodeMcadIdxs(allMcadIdxs);

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


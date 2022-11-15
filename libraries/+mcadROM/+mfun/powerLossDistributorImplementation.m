function nodeloss = powerLossDistributorImplementation(Tnodes, lossEachType, ...
    TrefStator, StatorCopperTempCoefResistivity, TrefRotor, RotorCopperTempCoefResistivity, LossDistrForEachType, numStates, CapMat, CoolantArrayIdxs) %#codegen
    % POWERLOSSDISTRIBUTORIMPLEMENTATION Calculates node loss based on the
    % node temperatures, loss of each type, and the loss distribution
    % table, amongst other parameters. It applies the temperature
    % correction for copper losses.
    % 
    % Input arguments:
    % - Tnodes: [double Nx1]: vector of node temperatures
    % - lossEachType: [double Mx1]: vector of loss value for each type of loss.
    % - TrefStator: [double 1x1]: Stator reference temperature.
    % - StatorCopperTempCoefResistivity: [double 1x1]: Stator temperature coefficient for copper resistivity
    % - TrefRotor: [double 1x1]: Rotor reference temperature.
    % - RotorCopperTempCoefResistivity: [double 1x1]: Rotor temperature coefficient for copper resistivity
    % - LossDistrForEachType [double MxN]: Distribution of node losses for each loss type
    % - numStates [double 1x1]: Number of states (nodes)
    % - CapMat [double Nx1]: Node thermal capacitance vector.
    % - CoolantArrayIdxs [double Px1]: Array of coolant indices.
    % 
    % Output arguments:
    % - nodeloss: [double Nx1]: Node power losses.

    % Copyright 2022 The MathWorks, Inc.
    
    % Stator copper loss temperature correction
    idxStatorCopperLoss = 1; 
    LossDistrForStatorCopperLoss = LossDistrForEachType(idxStatorCopperLoss,:); 
    statorCopperLoss = lossEachType(idxStatorCopperLoss);
    statorCopperLossNodes = LossDistrForStatorCopperLoss'*statorCopperLoss;
    statorCopperNodes = LossDistrForStatorCopperLoss(:)>0; 
    statorCopperNodes(CoolantArrayIdxs) = false; % Do not include coolant nodes
    TstatorAvg = sum(Tnodes(statorCopperNodes).*CapMat(statorCopperNodes)/sum(CapMat(statorCopperNodes))); 
    temperatureFactorStatorCopper = (1 + (TstatorAvg-20)*StatorCopperTempCoefResistivity)/(1 + (TrefStator-20)*StatorCopperTempCoefResistivity);

    % Rotor copper loss temperature correction
    idxRotorCopperLoss = 2; 
    LossDistrForRotorCopperLoss = LossDistrForEachType(idxRotorCopperLoss,:); 
    rotorCopperLoss = lossEachType(idxRotorCopperLoss);
    rotorCopperLossNodes = LossDistrForRotorCopperLoss'*rotorCopperLoss;
    rotorCopperNodes = LossDistrForRotorCopperLoss(:)>0; 
    rotorCopperNodes(CoolantArrayIdxs) = false; % Do not include coolant nodes
    TrotorAvg = sum(Tnodes(rotorCopperNodes).*CapMat(rotorCopperNodes)/sum(CapMat(rotorCopperNodes))); 
    temperatureFactorRotorCopper = (1 + (TrotorAvg-20)*RotorCopperTempCoefResistivity)/(1 + (TrefRotor-20)*RotorCopperTempCoefResistivity);
    
    % Non-temperature-dependent loss distributions
    otherlossNodes = zeros(numStates,1); 
    for idxLoss = 1:length(lossEachType) 
        if idxLoss ~= idxStatorCopperLoss && idxLoss ~= idxRotorCopperLoss % ignore copper loss types 
            otherlossNodes = otherlossNodes + LossDistrForEachType(idxLoss,:)'*lossEachType(idxLoss);
        end
    end
    
    nodeloss = zeros(numStates,1); 
    nodeloss(:) = statorCopperLossNodes.*temperatureFactorStatorCopper + rotorCopperLossNodes.*temperatureFactorRotorCopper + otherlossNodes; 

    nodeloss(1) = 0; % Ambient node is constant temperature, loss has no sense. 
    nodeloss(CoolantArrayIdxs) = 0; % Fluid heatflow to be computed in the State-Space model with upstream correction. 
end



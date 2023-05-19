classdef tPowerLossDistributorPassiveCooling < matlab.unittest.TestCase
    % Unit tests for PowerLosDistributor (Passive Cooling) library block

    % Copyright 2023 The MathWorks, Inc.

    properties
        modelName = 'mPowerLossDistributorPassiveCooling'
    end

    properties(TestParameter)
        TsimVal = {20, -10, 80}
        lossTypeVarName = {'armatureCopperLossVal', 'backIronLossVal', 'toothIronLossVal', 'magnetLossVal'};
    end
    
     methods(Test, ParameterCombination='pairwise')

         function testLossDistrDataConsitency(test)
            % Check that LossDistrForEachType is conservative i.e. sum of
            % rows = 1

             model = test.modelName;
             load_system(model);
             test.addTeardown(@close_system,model,0)

             mdlwks = get_param(model, 'ModelWorkspace');
             LossDistrForEachType = mdlwks.getVariable('LossDistrForEachType');
             Tgiven = mdlwks.getVariable('Tgiven');
             assignin(mdlwks, 'Tsim', Tgiven);

             [numLossTypes, ~] = size(LossDistrForEachType);

             for idxLossType = 1:numLossTypes
                thisLossDistr = LossDistrForEachType(idxLossType, :);
                if sum(thisLossDistr)>1e-3
                    test.verifyEqual(sum(thisLossDistr), 1, 'AbsTol', 1e-3, 'The loss distribution is not conservative');
                end
             end

         end

         function testTotalPowerConservation(test)
            % Check that the sum of input losses equals the sum of output
            % losses

             model = test.modelName;
             load_system(model);
             test.addTeardown(@close_system,model,0)

             mdlwks = get_param(model, 'ModelWorkspace');
             Tgiven = mdlwks.getVariable('Tgiven');
             assignin(mdlwks, 'Tsim', Tgiven);

             out = sim(model);

             lossEachType = out.yout{1}.Values.Data;
             nodeloss = out.yout{3}.Values.Data(:,1,end); % last timestamp

             test.verifyEqual(sum(lossEachType), sum(nodeloss), 'AbsTol', 1e-3, 'total losses in is not equal to total losses out');

         end

         function testIndividualLossTypes(test, lossTypeVarName)
            % For each loss type, check that loss in = loss out

             model = test.modelName;
             load_system(model);
             test.addTeardown(@close_system,model,0)

             mdlwks = get_param(model, 'ModelWorkspace');
             Tgiven = mdlwks.getVariable('Tgiven');
             assignin(mdlwks, 'Tsim', Tgiven);

             lossTypeVal = mdlwks.getVariable(lossTypeVarName);

             allLossTypesVarNames = test.lossTypeVarName;
             for idx = 1:length(allLossTypesVarNames)
                thisLossType = allLossTypesVarNames{idx};
                if ~strcmp(lossTypeVarName, thisLossType)
                    assignin(mdlwks, thisLossType, 0); % only lossTypeVarName is non-zero
                end
             end

             out = sim(bdroot);

             lossEachType = out.yout{1}.Values.Data;
             nodeloss = out.yout{3}.Values.Data(:,1,end); % last timestamp

             test.verifyEqual(sum(lossEachType), lossTypeVal, 'AbsTol', 1e-3, 'Loss each type does not sum to the value of the only non-zero loss');
             test.verifyEqual(sum(lossEachType), sum(nodeloss), 'AbsTol', 1e-3, 'total losses in is not equal to total losses out');

         end

         function testCopperLossTemperatureDependency(test, TsimVal)
            % Check that copper loss is correctly scaled based on the
            % copper resistivity temperature coefficient.

             model = test.modelName;
             load_system(model);
             test.addTeardown(@close_system,model,0)

             mdlwks = get_param(model, 'ModelWorkspace');
             assignin(mdlwks, 'Tsim', TsimVal);

             copperLossVal = mdlwks.getVariable('armatureCopperLossVal');

             allLossTypesVarNames = test.lossTypeVarName;
             for idx = 1:length(allLossTypesVarNames)
                thisLossType = allLossTypesVarNames{idx};
                if ~strcmp('armatureCopperLossVal', thisLossType)
                    assignin(mdlwks, thisLossType, 0); % only lossTypeVarName is non-zero
                end
             end

             out = sim(bdroot);

             lossEachType = out.yout{1}.Values.Data;
             nodeloss = out.yout{3}.Values.Data(:,1,end); % last timestamp

             Tgiven = mdlwks.getVariable('Tgiven');
             StatorCopperTempCoefResistivity = mdlwks.getVariable('StatorCopperTempCoefResistivity');

             temperatureFactorStatorCopper = (1 + (TsimVal-20)*StatorCopperTempCoefResistivity)/(1 + (Tgiven-20)*StatorCopperTempCoefResistivity);
             adjustedLossExpected = copperLossVal*temperatureFactorStatorCopper ;

             test.verifyEqual(sum(lossEachType), copperLossVal, 'AbsTol', 1e-3, 'Loss each type does not sum to the value of the input copper loss');
             test.verifyEqual(sum(nodeloss), adjustedLossExpected, 'AbsTol', 1e-3, 'nodeloss does not sum to the temperature-adjusted copper loss value');
             
         end

     end
end


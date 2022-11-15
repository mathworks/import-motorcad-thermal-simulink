classdef tThermalInterface < matlab.unittest.TestCase
    % Tests for mcadinterface.ThermalInterface

    % Copyright 2022 The MathWorks, Inc.
    
    properties
        realMotFile = {'data/e5_IM_HWJ.mot'}
        slkModel = {'mStateSpaceVsMcad', 'mLPVvsMcad'}
        objectUnderTest
        classUnderTest = @mcadinterface.ThermalInterface
    end
    
    methods(TestClassSetup)
        function classSetup(test)
            theProject = matlab.project.rootProject;
            projectRootDir = theProject.RootFolder;
            test.objectUnderTest = test.classUnderTest(fullfile(projectRootDir,test.realMotFile));
        end
    end

    % Unit tests for methods in mcadinterface.ThermalInterface
    methods(Test)

        function testCalculateThermalSteadyState(test)
            % Unit test for calculateThermalSteadyState method.
            % Checks that getting the ambient temperature node produces no
            % warning, which means a steady-state calculation has been
            % performed successfully.

            test.objectUnderTest.calculateThermalSteadyState();
            test.verifyWarningFree(@() test.objectUnderTest.getSteadyStateTemperatureForNodeMcadIdxs(0));
        end

        function testCalculateThermalTransient(test)
            % Unit test for calculateThermalTransient method.
            % Checks that getting the ambient node transient temperature produces no
            % warning, which means a transient calculation has been
            % performed successfully.

            test.objectUnderTest.calculateThermalTransient();
            test.verifyWarningFree(@() test.objectUnderTest.getTransientTemperatureForNodeMcadIdxs(0));
        end

        function testRunThermalSteadyStateWithSpecifiedLosses(test)
            % Unit test for runThermalSteadyStateWithSpecifiedLosses method
            % Checks that the TnodesVec output has the expected size.

            lossVec = transpose(1:16); % [W]
            test.objectUnderTest.runThermalSteadyStateWithSpecifiedLosses(lossVec);

            % Consistency check
            actualLossValues = test.objectUnderTest.LossValues;
            test.verifyEqual(actualLossValues, lossVec, 'AbsTol', 1e-3, 'LossValues does not match expected values');

            allMcadIdxs = [test.objectUnderTest.NodeNamesAndMcadIdx{:,2}];           
            TnodesVec =  test.objectUnderTest.getSteadyStateTemperatureForNodeMcadIdxs(allMcadIdxs);

            test.verifyEqual(length(TnodesVec), length(allMcadIdxs), ...
                        'TnodesVec does not have the expected size');
        end

        function testRunThermalSteadyStateWithSpecifiedTorqueSpeed(test)
            % Unit test for runThermalSteadyStateWithSpecifiedTorqueSpeed method
            % Checks that the TnodesVec output has the expected size.

            torqueVal = 50; % N*m
            speedVal = 3000; % rpm
            test.objectUnderTest.runThermalSteadyStateWithSpecifiedTorqueSpeed(torqueVal, speedVal);

            allMcadIdxs = [test.objectUnderTest.NodeNamesAndMcadIdx{:,2}];           
            TnodesVec =  test.objectUnderTest.getSteadyStateTemperatureForNodeMcadIdxs(allMcadIdxs);

            test.verifyEqual(length(TnodesVec), length(allMcadIdxs), ...
                        'TnodesVec does not have the expected size');
        end

        function testRunThermalTransientWithSpecifiedLosses(test)
            % Unit test for runThermalTransientWithSpecifiedLosses method
            % Checks that the TnodesVecMcad output has the expected size.

            lossVec = 30*ones(16,1);
            test.objectUnderTest.EnableStatorTempCoeffRes = 1;
            test.objectUnderTest.EnableRotorTempCoeffRes = 1;
            stopTime = 500;
            numTimeSteps = 50;
            test.objectUnderTest.runThermalTransientWithSpecifiedLosses(lossVec, stopTime, numTimeSteps);

            allMcadIdxs = [test.objectUnderTest.NodeNamesAndMcadIdx{:,2}];
            [tVec, TnodesVecMcad] = test.objectUnderTest.getTransientTemperatureForNodeMcadIdxs(allMcadIdxs);

            test.verifyEqual(tVec(end), stopTime, 'Final time step is not equal to stopTime');
            test.verifyEqual(size(TnodesVecMcad), [length(allMcadIdxs), length(tVec)], ...
                        'TnodesVecMcad does not have the expected size');

        end

        function testRunThermalTransientWithSpecifiedTorqueSpeed(test)
            % Unit test for runThermalTransientWithSpecifiedTorqueSpeed method
            % Checks that the TnodesVecMcad output has the expected size.

            torqueVal = 50; % N*m
            speedVal = 3000; % rpm
            stopTime = 500;
            numTimeSteps = 50;
            test.objectUnderTest.runThermalTransientWithSpecifiedTorqueSpeed(torqueVal, speedVal, stopTime, numTimeSteps);

            allMcadIdxs = [test.objectUnderTest.NodeNamesAndMcadIdx{:,2}];           
            [tVec, TnodesVecMcad] = test.objectUnderTest.getTransientTemperatureForNodeMcadIdxs(allMcadIdxs);

            
            test.verifyEqual(tVec(end), stopTime, 'Final time step is not equal to stopTime');
            test.verifyEqual(size(TnodesVecMcad), [length(allMcadIdxs), length(tVec)], ...
                        'TnodesVecMcad does not have the expected size');
        end
      
        function testUpdateMatricesAndGroupNamesAndNodes(test)
            % Unit test for updateMatricesAndGroupNamesAndNodes method
            % Checks that the thermal matrices are correctly updated.

            test.objectUnderTest.updateMatricesAndGroupNamesAndNodes();

            theProject = matlab.project.rootProject;
            projectRootDir = theProject.RootFolder;
            motFullFile = fullfile(projectRootDir,test.realMotFile);
            [CapMat, ResMat, PowMat, TempMat, ~] = getThermalMatricesFromFiles(motFullFile);
            % Basic checks
            relTol = 2e-3;
            test.verifyEqual(CapMat, test.objectUnderTest.CapMat,'RelTol', relTol);
            test.verifyEqual(ResMat, test.objectUnderTest.ResMat,'RelTol', relTol);
            test.verifyEqual(PowMat, test.objectUnderTest.PowMat,'RelTol', relTol);
            test.verifyEqual(TempMat, test.objectUnderTest.TempMat,'RelTol', relTol);

        end

        function testGetSteadyStateTemperatureForNodeMcadIdxs(test)
            % Unit test for getSteadyStateTemperatureForNodeMcadIdxs method
            % Checks that the ambient temperature node has the expected

            AmbientNodeMcadIdx = 0;
            test.objectUnderTest.Tambient_degC = 37;
            test.objectUnderTest.calculateThermalSteadyState();
            Tambient = test.objectUnderTest.getSteadyStateTemperatureForNodeMcadIdxs(AmbientNodeMcadIdx);
            test.verifyEqual(Tambient, test.objectUnderTest.Tambient_degC);
        end

        function testGetSteadyStatePowerForNodeMcadIdxs(test)
            % Unit test for getSteadyStatePowerForNodeMcadIdxs method
            % Checks that the ambient node power is zero, as expected.

            AmbientNodeMcadIdx = 0;
            test.objectUnderTest.calculateThermalSteadyState();
            Pambient = test.objectUnderTest.getSteadyStatePowerForNodeMcadIdxs(AmbientNodeMcadIdx);
            test.verifyEqual(Pambient, 0); % Mcad convention: fixed temperature has 0 power
        end

        function testUpdateCoolingSystemsData(test)
            % Unit test for updateCoolingSystemsData method
            % Checks that the cooling system data is correct when one
            % cooling system is enabled.

            test.objectUnderTest.updateCoolingSystemsData();
            test.verifyEqual(test.objectUnderTest.CoolingSystemNamesAndMcadIdxes, ...
                    {'Housing Water Jacket', [29, 132]}, 'Actual cooling system cell is not equal to expected value');

        end

        function testUpdateCoolingSystemDataWithTwoCoolingSys(test)
            % Unit test for updateCoolingSystemsData method
            % Checks that the cooling system data is correct when two
            % cooling systems are enabled.

            test.objectUnderTest.Ventilated_Enable = 1;
            test.objectUnderTest.updateMatricesAndGroupNamesAndNodes();
            test.objectUnderTest.updateCoolingSystemsData();

            test.verifyEqual(test.objectUnderTest.CoolingSystemNamesAndMcadIdxes, ...
                    {'Ventilated', [52 53 54 55 56 57 59 60 98 99 131]; ...
                     'Housing Water Jacket', [29 132] ...
                     }, 'Actual cooling system cell is not equal to expected value');

        end

        function testGetAdjacencyMatAndInletOutletIdxs(test)
            % Unit test for getAdjacencyMatAndInletOutletIdxs method
            % Checks that the outputs of getAdjacencyMatAndInletOutletIdxs
            % match the expected values.

            test.objectUnderTest.Ventilated_Enable = 0;
            test.objectUnderTest.HousingWaterJacket_Enable = 1;
            test.objectUnderTest.ShaftSpiralGroove_Enable = 0;
            % Cooling inlet values
            test.objectUnderTest.HousingWaterJacket_FlowRate_m3ps = 5/1000/60; % 5 lpm
            test.objectUnderTest.HousingWaterJacket_InletTemperature_degC = 30;
            % Update
            test.objectUnderTest.updateModel();

            [AdjacencyMat, InletArrayIdxs, OutletArrayIdxs, CoolantArrayIdxs] = test.objectUnderTest.getAdjacencyMatAndInletOutletIdxs();

            test.verifyEqual(InletArrayIdxs, 26, 'Inlet indexes do not match expected value');
            test.verifyEqual(OutletArrayIdxs, 21, 'Outlet indexes do not match expected value');
            test.verifyEqual(CoolantArrayIdxs, [21; 26], 'Coolant indexes do not match expected value');
            test.verifyEqual(AdjacencyMat(26,21), 1, 'Connectivity of adjacency mat for nodes 29->132 is not correct');
            test.verifyEqual(AdjacencyMat(21,26), 0, 'Connectivity of adjacency mat for nodes 29->132 is not correct');

        end

        function testGenerateSimulinkReducedOrderModel(test)
            % Test for generateSimulinkReducedOrderModel method.
            % Checks that the generated model data is consistent, and checks
            % that the maximum temperature error in a transient simulation
            % compared with Motor-CAD is acceptably low for two different table
            % breakpoints.
            
            modelName = 'tempGeneratedModel';
            bdclose(modelName); % in case it exists
            modelNamePath = findModelInProject(modelName);
            if ~isempty(modelNamePath)
                delete(modelNamePath);
            end
            coolingSystemsEnabled = {'Housing Water Jacket'};
            BkptsStruct = struct();
            BkptsStruct.w = [1000,3000]; % rpm
            BkptsStruct.fr1 = [5,10]; % lpm
            BkptsStruct.Tin1 = [20,30]; % degC
            
            [AmatND, BmatND, ResMatND, TnodesInit, AdjacencyMat, ~, ~] = ...
                    test.objectUnderTest.generateSimulinkReducedOrderModel(modelName, coolingSystemsEnabled, BkptsStruct);

            test.verifyEqual(size(AmatND), [2,2,2,150,150], 'Size of A state-space matrix not equal to expected value');
            test.verifyEqual(size(BmatND), [2,2,2,150,150], 'Size of B state-space matrix not equal to expected value');
            test.verifyEqual(size(ResMatND), [2,2,2,150,150],'Size of resistance matrix not equal to expected value');
            test.verifyEqual(size(TnodesInit), [150,1],'Size of TnodesInit not equal to expected value');
            test.verifyEqual(size(AdjacencyMat), [150,150],'Size of AdjacencyMat not equal to expected value');
            
            torVal = 50; % N*m
            maxTempError111 = compareSimulinkAndMcadAtBkpt(modelName, test.realMotFile{1}, [1,1,1], torVal);
            maxTempError222 = compareSimulinkAndMcadAtBkpt(modelName, test.realMotFile{1}, [2,2,2], torVal);
            Ttol = 2; % degC
            test.verifyLessThan(maxTempError111, Ttol, 'Max temperature error is not acceptable for [1,1,1] breakpoints');
            test.verifyLessThan(maxTempError222, Ttol, 'Max temperature error is not acceptable for [2,2,2] breakpoints');

            bdclose(modelName);

        end
    end
end
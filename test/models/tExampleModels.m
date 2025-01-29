classdef tExampleModels < matlab.unittest.TestCase
    % Verify the example models included in this repository run with no
    % warnings or errors

    % Copyright 2025 The MathWorks, Inc.

    properties(TestParameter)
        modelsToTest = {
            'e5_IM_HWJ_ROM'; 
            'e8_IPMSM_HWJandVent_ROM'; 
            'DriveCycle_e8_IPMSM_HWJandVent_ROM'; 
            'SimscapeSystemLevel_e8_IPMSM_HWJandVent_ROM';
        }
    end

    methods(Test)
        function checkRunWarnFree(testCase, modelsToTest)
            % Get the exact path to load the model
            modelNamePath = findModelInProject(modelsToTest);
            testCase.verifyNotEmpty(modelNamePath, ...
                sprintf('Model "%s" not found in the project.', modelsToTest));

            % Load the model
            load_system(modelNamePath);
            testCase.addTeardown(@close_system,modelNamePath,0)

            % Verify it runs warning-free (just one second of simulation)
            testCase.verifyWarningFree(@() sim(modelsToTest, "StopTime", "1"), ...
                sprintf('Model "%s" generated warnings during simulation.', modelsToTest));
        end
    end
end
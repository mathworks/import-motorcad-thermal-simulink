classdef tLossMapsWithDCV < matlab.unittest.TestCase
    % Unit tests for one 3-D Lookup Table version of LossMapsWithDCV block
    % (Torque × Speed × DC-Voltage → Losses)
    %
    % Copyright 2025 The MathWorks, Inc.

    %% Test-fixture data --------------------------------------------------
    properties
        % Name of the Simulink model under test. The model is expected to
        % take the vectors xVec (torque), yVec (speed), zVec (DC-voltage)
        % and zg3D (loss grid) from the base workspace as well as the
        % scalars torqueVal, speedVal and voltVal that define the lookup
        % point used in this test.
        modelName = "mOneLossMap3D";
    end

    %% Simple 3-D cases ---------------------------------------------------
    % Each CaseData structure contains vectors that define the three axes
    % and a 3-D array "LossArray" whose dimensions correspond to
    %   size(TorqueVec) × size(SpeedVec) × size(VoltageVec)
    %
    % For readability only one representative case is supplied here.
    properties (TestParameter)
        CaseData = { ...
            struct( ...
                "ShaftTorqueVec",  [-1; 0; 1], ...               %   3 × 1
                "SpeedVec",        [0 1e3 2e3], ...              %   1 × 3
                "VoltageVec",      [200 400], ...                %   1 × 2
                "LossArray",       cat( 3, ...
                    200 * 1e-3 * repmat([0 1e3 2e3], 3, 1).^2 .* repmat([-1; 0; 1], 1, 3), ...
                    400 * 1e-3 * repmat([0 1e3 2e3], 3, 1).^2 .* repmat([-1; 0; 1], 1, 3)  ...
                ) ...
            ) ...
        };
    end

    %% Test methods -------------------------------------------------------
    methods (Test, ParameterCombination = "pairwise")
        function test3DInterpolation(test, CaseData)
            % Verify that the 3-D Lookup Table block reproduces the map
            % values at grid points (no extrapolation involved).

            % Bring Simulink model into memory and ensure it is closed
            % afterwards even if the test fails.
            mdl = test.modelName;
            load_system(mdl);
            test.addTeardown(@close_system, mdl, 0);

            % Extract the case-specific data --------------------------------
            TorqueVec   = CaseData.ShaftTorqueVec;
            SpeedVec    = CaseData.SpeedVec;
            VoltageVec  = CaseData.VoltageVec;
            LossArray   = CaseData.LossArray;

            % Prepare lookup table vectors and grid for the model ----------
            % Convert the axis vectors to a full 3-D mesh 
            [TorqueMat, SpeedMat, VoltageMat] = ndgrid( ...
                TorqueVec, SpeedVec, VoltageVec);
            % reinterpolate
            [xVec, yVec, zVec, zg3D] = reInterpolateTable3D( ...
                TorqueMat, SpeedMat, VoltageMat, LossArray);


            assignin("base", "xVec",    xVec);    % torque axis
            assignin("base", "yVec",    yVec);    % speed  axis
            assignin("base", "zVec",    zVec);    % voltage axis

            % Publish grid under both possible names expected by models ----
            assignin("base", "zg3D",    zg3D);    % new tests
            assignin("base", "zgMat",   zg3D);    % legacy model alias

            % Choose the corner value (last indices) to check --------------
            torqueVal   = TorqueVec(end);
            speedVal    = SpeedVec(end);
            voltVal     = VoltageVec(end);

            assignin("base", "torqueVal",  torqueVal);
            assignin("base", "speedVal",   speedVal);

            % Publish lookup point under both spellings so all blocks work --
            assignin("base", "voltVal",    voltVal);   % used by this test
            assignin("base", "voltageVal", voltVal);   % used by model

            % Run the model -------------------------------------------------
            out = sim(mdl);

            % Retrieve the last sample of the output signal -----------------
            zActual   = out.yout{1}.Values.Data(end);
            zExpected = LossArray(end, end, end);

            % Verification --------------------------------------------------
            test.verifyEqual(zActual,   zExpected, ...
                "Actual interpolated value does not match expected value");
            test.verifyEqual(zg3D,      LossArray, ...
                "Gridded loss table does not match the expected map");
        end
    end
end

%% Helper function ---------------------------------------------------------
function [xVec, yVec, zVec, vGrid] = reInterpolateTable3D(xMat, yMat, zMat, vMat)
%REINTERPOLATETABLE3D  Re-interpolate scattered data onto a cubic grid
%                      (same policy as the 2-D helper `reInterpolateTable`).
%
%   INPUTS
%     xMat, yMat, zMat  – equal-sized 3-D arrays of sample coordinates
%     vMat              – same-sized 3-D array of values at those points
%
%   OUTPUTS
%     xVec, yVec, zVec  – evenly-spaced axis vectors
%     vGrid             – |numel(xVec)|-by-|numel(yVec)|-by-|numel(zVec)| grid

    F = scatteredInterpolant( ...
            xMat(:), yMat(:), zMat(:), vMat(:), ...
            "linear", "nearest");

    [nx, ny, nz] = size(xMat);          % identical for all four inputs
    xVec = linspace(min(xMat(:)), max(xMat(:)), nx);
    yVec = linspace(min(yMat(:)), max(yMat(:)), ny);
    zVec = linspace(min(zMat(:)), max(zMat(:)), nz);

    [X, Y, Z] = ndgrid(xVec, yVec, zVec);
    vGrid     = F(X, Y, Z);
end
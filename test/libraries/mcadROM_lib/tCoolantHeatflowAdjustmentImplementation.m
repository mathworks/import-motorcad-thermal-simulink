classdef tCoolantHeatflowAdjustmentImplementation < matlab.unittest.TestCase
    % Unit tests for mcadROM.mfun.coolantHeatflowAdjustmentImplementation

    % Copyright 2022 The MathWorks, Inc.

    properties
        funcUnderTest = @mcadROM.mfun.coolantHeatflowAdjustmentImplementation
    end

     methods(Test)

         function testCoolantHeatflowAdjustmentSimpleCase(test)
            % Simple case with 4 nodes and 3 coolant nodes

            numStates = 4;
            Tnodes = [100,200,300,400]';
            ResMat = 1e9*ones(4,4);
            for idx = 1:4
                ResMat(idx,idx) = 0;
            end
            % Coolant connection path: 1->2->3
            CoolAdjMat = zeros(4,4);
            CoolAdjMat(1,2) = 1;
            CoolAdjMat(2,3) = 1;
            CoolantArrayIdxs = [1,2,3];
            ResMat(1,2) = 1; ResMat(2,1) = 1;
            ResMat(2,3) = 1; ResMat(3,2) = 1;

            nodelossCorrectionActual = test.funcUnderTest(Tnodes, ResMat, CoolAdjMat, CoolantArrayIdxs, numStates);

            % Compute expected nodelossCorrection vector
            nodelossCorrectionExpected = zeros(size(Tnodes));
            nodelossCorrectionExpected(1) = -(Tnodes(2) - Tnodes(1))/ResMat(1,2);
            nodelossCorrectionExpected(2) = -(Tnodes(3) - Tnodes(2))/ResMat(2,3);
            nodelossCorrectionExpected(3) = 0;

            test.verifyEqual(nodelossCorrectionActual, nodelossCorrectionExpected, ...
                                'Node loss correction does not match expected values');
            
         end

     end
end


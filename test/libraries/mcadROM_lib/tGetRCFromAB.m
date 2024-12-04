classdef tGetRCFromAB < matlab.unittest.TestCase
    % Unit tests for mcadROM.mfun.getRCFromAB
    
    % Copyright 2022 The MathWorks, Inc.
    
    properties
        funPreProcess = @getStateSpaceMatricesFromThermalMatrices
        funcUnderTest = @mcadROM.mfun.getRCFromAB
    end

     methods(Test)

         function testConsistencyRC2RC(test)
            
             % Check that the function returns the original data,
             % generated from getStateSpaceMatricesFromThermalMatrices
            
             CapMat1 = 1e-2*[1,2,3,4]';
             ResMat1 = rand(4);
             for idx=1:4
                 ResMat1(idx,idx) = 0;
             end
             ResMat1=triu(ResMat1.',1) + tril(ResMat1); % make symmetrical
            
             test.assertEqual(ResMat1, ResMat1', 'Resistance matrix is not symmetrical');

             [Amat, Bmat] = test.funPreProcess(CapMat1, ResMat1);
             [ResMat2, CapMat2] = test.funcUnderTest(Amat,Bmat);

             test.verifyEqual(ResMat1, ResMat2, 'AbsTol', 1e-7, 'Resistance matrix not equal to the expected matrix');
             test.verifyEqual(CapMat1, CapMat2, 'AbsTol', 1e-7, 'Capacitance matrix not equal to the expected matrix');
             
         end

          function testConsistencyAB2AB(test)
            
             % Check that the function returns the original data,
             % generated from getStateSpaceMatricesFromThermalMatrices
            
             Amat1 = ones(3);
             Amat1(1,1) = 20;
             Amat1(2,2) = 30;
             Amat1(3,3) = 40;
             Bmat1 = diag(rand(3,1));
            
             [ResMat1, CapMat1] = test.funcUnderTest(Amat1,Bmat1);
             [Amat2, Bmat2] = test.funPreProcess(CapMat1, ResMat1);
             [ResMat2, CapMat2] = test.funcUnderTest(Amat2,Bmat2);

             test.verifyEqual(ResMat1, ResMat2, 'AbsTol', 1e-7, 'Resistance matrix not equal to the expected matrix');
             test.verifyEqual(CapMat1, CapMat2, 'AbsTol', 1e-7, 'Capacitance matrix not equal to the expected matrix');
             

         end
     end
end

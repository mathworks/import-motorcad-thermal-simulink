classdef tLossMaps < matlab.unittest.TestCase
    % Unit tests for 2D Lookup Table block in the LossMaps library block

    % Copyright 2022 The MathWorks, Inc.
    
    properties
        modelName = 'mOneLossMap';
    end

    % Simple cases
    properties(TestParameter)
        CaseData = { struct('ShaftTorqueMat', repmat([-1,0,1]',1,3), 'SpeedMat', repmat([0,1000,2000],3,1), 'LossMat', 1e-3*repmat([0,1000,2000],3,1).^2.*repmat([-1,0,1]',1,3)), ...
                     struct('ShaftTorqueMat', repmat([-1,0,1,2]',1,3), 'SpeedMat', repmat([0,1000,2000],4,1), 'LossMat', 1e-3*repmat([0,1000,2000],4,1).^2.*repmat([-1,0,1,2]',1,3)), ...
                     }
    end
    
     methods(Test, ParameterCombination='pairwise')

         function testInterpolation(test, CaseData)
             % Check that the 2D Lookup Table blocks correctly interpolate
             % the maps
            
             model = test.modelName;
             load_system(model);
             test.addTeardown(@close_system,model,0)

             ShaftTorqueMat = CaseData.ShaftTorqueMat;
             SpeedMat = CaseData.SpeedMat;
             LossMat = CaseData.LossMat;

             [xVec, yVec, zgMat] = reInterpolateTable(ShaftTorqueMat, SpeedMat, LossMat);
             assignin('base', 'xVec', xVec);
             assignin('base', 'yVec', yVec);
             assignin('base', 'zgMat', zgMat);

             torqueVal = ShaftTorqueMat(end,end);
             speedVal = SpeedMat(end,end);
             assignin('base', 'torqueVal', torqueVal);
             assignin('base', 'speedVal', speedVal);

             out = sim(model);

             zActual = out.yout{1}.Values.Data(end);
             zExpected = LossMat(end,end);

             test.verifyEqual(zActual,zExpected, 'Actual interpolated value does not match expected value');
             test.verifyEqual(zgMat, LossMat, 'Actual table does not match expected value');

         end
     
     end
end



function [xVec, yVec, zgMat] = reInterpolateTable(xsMat, ysMat, zsMat)
        % Returns grid-interpolated table from scattered table.

        zInterpolant = scatteredInterpolant(xsMat(:), ysMat(:), zsMat(:), 'linear', 'nearest');
        xMin = min(xsMat(:));
        xMax = max(xsMat(:));
        [xLen,~] = size(xsMat);
        xVec = linspace(xMin, xMax, xLen);
        yMin = min(ysMat(:));
        yMax = max(ysMat(:));
        [~,yLen] = size(ysMat);
        yVec = linspace(yMin, yMax, yLen);
        [xgrid, ygrid] = ndgrid(xVec, yVec);
        zgMat = zInterpolant(xgrid, ygrid);

end
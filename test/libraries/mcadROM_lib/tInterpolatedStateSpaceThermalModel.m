classdef tInterpolatedStateSpaceThermalModel < matlab.unittest.TestCase
    % tInterpolatedStateSpaceThermalModel Tests the Simulink subsystem
    % InterpolatedStateSpaceThermalModel
    % This also indirectly tests the function getStateSpaceMatricesFromThermalMatrices

    % Copyright 2022 The MathWorks, Inc.
    
    properties
        modelName = {'mSimpleThermalStateSpaceLPV', 'mLPVsys'};
        funcPreProcess = @getStateSpaceMatricesFromThermalMatrices
    end

    % Simple cases
    properties(TestParameter)
        CaseData = { struct('CapMat', [1; 1; 1], 'ResMat',  [1,1,1; 1,1,1; 1,1,1]), ...
                     struct('CapMat', [1e-1; 10; 1e2], 'ResMat',  [1,1e3,1; 1e3,1,1; 1,1,1]), ...
                     struct('CapMat', [1; 1e4; 1], 'ResMat',  [1,1,1e4; 1,1,1; 1e4,1,1])}

        PowData = {[0;0;0], ... % [W]
                   [1000;500;2000], ...
                   [500;1e4;0],...
                   [0;0;-100]}

        dimInterp = {1,2}
        wBkpt = {1,2}
        frBkpt = {1,2}

        idxSScomb = {[1,1,1],[2,2,2],[1,2,1],[2,1,2]}
        nodelossVal = {2*ones(150,1), 4*rand(150,1)}
    end
    
     methods(Test, ParameterCombination='exhaustive')

         function testLPVvsStateSpace(test, idxSScomb, nodelossVal)
            % Verify A,B matrices and node temperatures from the LPV block 
            % match the values from the State Space block, for different
            % index combinations and nodeloss values.

             import Simulink.sdi.constraints.MatchesSignal

             model = test.modelName{2};
             load_system(model);
             test.addTeardown(@close_system,model,0)

             mdlwks = get_param(model, 'ModelWorkspace');
             StateSpaceND = mdlwks.getVariable('StateSpaceND');
             samplingGrid = StateSpaceND.SamplingGrid;
             wBkpts = squeeze(samplingGrid.w(:,1,1))';
             frBkpts = squeeze(samplingGrid.fr1(1,:,1));
             TinBkpts = squeeze(samplingGrid.Tin1(1,1,:))';

             wIdx = idxSScomb(1);
             frIdx = idxSScomb(2);
             TinIdx = idxSScomb(3);
             ssIdxComb = StateSpaceND(:,:,wIdx,frIdx,TinIdx);
             Amat = ssIdxComb.A;
             Bmat = ssIdxComb.B;
             [numStates, ~] = size(Amat);
             TnodesInit = TinBkpts(TinIdx).*ones(numStates,1);

             clear(mdlwks); % to ensure it uses data generated in this test 
             assignin(mdlwks, 'StateSpaceND', StateSpaceND);
             assignin(mdlwks, 'Amat', Amat);
             assignin(mdlwks, 'Bmat', Bmat);
             assignin(mdlwks, 'speedVal', wBkpts(wIdx));
             assignin(mdlwks, 'flowrateVal', frBkpts(frIdx));
             assignin(mdlwks, 'TinVal', TinBkpts(TinIdx));
             assignin(mdlwks, 'TnodesInit', TnodesInit);
             assignin(mdlwks, 'nodelossVal', nodelossVal);

             out = sim(model);

             % A,B matrices verification
             AmatLPV = out.yout{1}.Values.Data(:,:,end);
             BmatLPV = out.yout{2}.Values.Data(:,:,end);
             test.verifyEqual(AmatLPV, Amat, 'AbsTol', 1e-8, 'RelTol', 1e-8, 'A matrix does not match expected value');
             test.verifyEqual(BmatLPV, Bmat, 'AbsTol', 1e-8, 'RelTol', 1e-8, 'B matrix does not match expected value');

             % Output verification (Temperatures)
             TnodesVecLPV = squeeze(out.yout{3}.Values.Data(:,1,:));
             TnodesLPV = timeseries(TnodesVecLPV', out.tout);
             TnodesVecSS = squeeze(out.yout{4}.Values.Data(:,1,:));
             TnodesSS = timeseries(TnodesVecSS', out.tout);

             TabsTol = 1e-1; % degC
             relTol = 1e-3;
             timeTol = 1e-1;
             test.verifyThat(TnodesLPV, MatchesSignal(TnodesSS, 'AbsTol', TabsTol, 'RelTol', relTol, 'TimeTol', timeTol), 'LPV Temperatures do not match SS Temperatures');     
             
         end

         function testSimple3dofCaseAgainstSimscapeBaseline(test, CaseData, PowData)
             % Check that the StateSpaceThermalModelSubsystem produces the
             % same results as the equivalent Simscape thermal model (using
             % thermal resistance and thermal capacitance blocks) for a
             % thermal network with 3 nodes, with constant A,B matrices

             CapMat = CaseData.CapMat;
             ResMat = CaseData.ResMat;
             test.assertEqual(ResMat(1,2),ResMat(2,1), 'Invalid data: Non-symmetric resistance matrix');
             test.assertEqual(ResMat(1,3),ResMat(3,1), 'Invalid data: Non-symmetric resistance matrix');
             test.assertEqual(ResMat(2,3),ResMat(3,2), 'Invalid data: Non-symmetric resistance matrix');
             test.assertGreaterThan(CapMat, 0, 'Invalid data: non-positive capacitance vector');
             assignin('base', 'CapMat', CapMat);
             assignin('base', 'ResMat', ResMat);

             [AMat, BMat] = test.funcPreProcess(CapMat, ResMat);
           
             model = test.modelName{1};
             load_system(model);
             test.addTeardown(@close_system,model,0)

             TcoolantIn = 300;
             assignin('base', 'TcoolantIn', TcoolantIn);
             BkptsStruct = struct();
             BkptsStruct.w = [0,2];
             BkptsStruct.fr1 = [0,2];
             BkptsStruct.Tin1 = [TcoolantIn-100,TcoolantIn+100];
             assignin('base', 'BkptsStruct', BkptsStruct);
             
             % Simulate a point in the middle of the grid [0,2]
             wVal = 1;
             frVal = 1;
             assignin('base', 'wVal', wVal);
             assignin('base', 'frVal', frVal);

             AMatND = nan(2,2,2,3,3);
             BMatND = nan(2,2,2,3,3);
             ResMatND = nan(2,2,2,3,3);
             for idx1 = 1:2
                 for idx2 = 1:2
                     for idx3 = 1:2
                         AMatND(idx1,idx2,idx3, :,:) = AMat;
                         BMatND(idx1,idx2,idx3, :,:) = BMat;
                         ResMatND(idx1,idx2,idx3, :,:) = ResMat;
                         StateSpaceND(:,:,idx1,idx2,idx3) = ss(AMat, BMat, eye(3), zeros(3)); %#ok<AGROW> 
                     end
                 end
             end
             samplingGridStruct = getSamplingGridFromBkpts(BkptsStruct);
             StateSpaceND.SamplingGrid = samplingGridStruct;

             test.assertFalse(isnan(any(AMatND(:))));
             test.assertFalse(isnan(any(BMatND(:))));
             assignin('base', 'AMatND', AMatND);
             assignin('base', 'BMatND', BMatND);
             assignin('base', 'ResMatND', ResMatND);
             assignin('base', 'StateSpaceND', StateSpaceND);
             assignin('base', 'ResMat', ResMat);
             TnodesInit = [TcoolantIn;330;290];
             assignin('base', 'TnodesInit', TnodesInit);
             AdjMat = zeros(3,3);
             PowVec = PowData; % W
             idx_TcoolantIn = 1;
             idx_TcoolantOut = 2;
             AdjMat(idx_TcoolantIn, idx_TcoolantOut) = 1;
             assignin('base','PowVec', PowVec);
             assignin('base','idx_TcoolantIn', idx_TcoolantIn);
             assignin('base','idx_TcoolantOut', idx_TcoolantOut);
             assignin('base', 'AdjMat', AdjMat);

             out = sim(model);
             
             TnodesExpected = squeeze(out.yout{1}.Values.Data);
             TnodesGridded = squeeze(out.yout{2}.Values.Data)';
             TnodesLPV = squeeze(out.yout{3}.Values.Data)';

             test.verifyEqual(TnodesGridded, TnodesExpected, 'AbsTol', 1e-4, 'RelTol', 1e-4, ...
                 'griddedInterpolant version of InterpolatedStateSpaceThermalModel does not match Simscape baseline');
             test.verifyEqual(TnodesLPV, TnodesExpected, 'AbsTol', 1e-4, 'RelTol', 1e-4, ...
                 'LPV version of InterpolatedStateSpaceThermalModel does not match Simscape baseline');

         end

         function testInterpolation(test, dimInterp, wBkpt, frBkpt)
            
             % Check that the StateSpaceThermalModelSubsystem produces the
             % same results as the equivalent Simscape thermal model (using
             % thermal resistance and thermal capacitance blocks) for a
             % thermal network with 3 nodes, with interpolated A,B
             % matrices.

             CapMat1 = [1; 1; 1];
             ResMat1 = [1,1,1; 1,1,1; 1,1,1];
             CapMat2 = [2; 1; 3];
             ResMat2 = 0.5*[1,1,1; 1,1,1; 1,1,1];

             [AMat1, BMat1] = test.funcPreProcess(CapMat1, ResMat1);
             [AMat2, BMat2] = test.funcPreProcess(CapMat2, ResMat2);
           
             model = test.modelName{1};
             load_system(model);
             test.addTeardown(@close_system,model,0)

             TcoolantIn = 300;
             assignin('base', 'TcoolantIn', TcoolantIn);
             BkptsStruct = struct();
             BkptsStruct.w = [0,2];
             BkptsStruct.fr1 = [0,2];
             BkptsStruct.Tin1 = [TcoolantIn-100,TcoolantIn+100];
             assignin('base', 'BkptsStruct', BkptsStruct);
            
             AMatND = nan(2,2,2,3,3);
             BMatND = nan(2,2,2,3,3);
             ResMatND = nan(2,2,2,3,3);
             for idx1 = 1:2
                 for idx2 = 1:2
                     for idx3 = 1:2
                         switch dimInterp
                             case 1
                                flagInterp = idx1;
                             case 2
                                flagInterp = idx2;
                         end
                         if flagInterp ==1 
                             % 1st bkpt -> AMat1, BMat1
                             AMatND(idx1,idx2,idx3, :,:) = AMat1;
                             BMatND(idx1,idx2,idx3, :,:) = BMat1;
                             ResMatND(idx1,idx2,idx3, :,:) = ResMat1;
                             StateSpaceND(:,:,idx1,idx2,idx3) = ss(AMat1, BMat1, eye(3), zeros(3)); %#ok<AGROW> 
                         else
                             % 2nd bkpt -> AMat2, BMat2
                             AMatND(idx1,idx2,idx3, :,:) = AMat2;
                             BMatND(idx1,idx2,idx3, :,:) = BMat2;
                             ResMatND(idx1,idx2,idx3, :,:) = ResMat2;
                             StateSpaceND(:,:,idx1,idx2,idx3) = ss(AMat2, BMat2, eye(3), zeros(3)); %#ok<AGROW> 
                         end
                     end
                 end
             end
             samplingGridStruct = getSamplingGridFromBkpts(BkptsStruct);
             StateSpaceND.SamplingGrid = samplingGridStruct;

             test.assertFalse(isnan(any(AMatND(:))));
             test.assertFalse(isnan(any(BMatND(:))));
             assignin('base', 'AMatND', AMatND);
             assignin('base', 'BMatND', BMatND);
             assignin('base', 'ResMatND', ResMatND);
             assignin('base', 'StateSpaceND', StateSpaceND);

             TnodesInit = [TcoolantIn;330;290];
             assignin('base', 'TnodesInit', TnodesInit);
             AdjMat = zeros(3,3);
             PowVec = [200;400;100]; % W
             idx_TcoolantIn = 1;
             idx_TcoolantOut = 2;
             AdjMat(idx_TcoolantIn, idx_TcoolantOut) = 1;
             assignin('base','PowVec', PowVec);
             assignin('base','idx_TcoolantIn', idx_TcoolantIn);
             assignin('base','idx_TcoolantOut', idx_TcoolantOut);
             assignin('base', 'AdjMat', AdjMat);

             % CASE 1: 1st bkpt -----
             assignin('base', 'CapMat', CapMat1);
             assignin('base', 'ResMat', ResMat1);
             switch dimInterp
                 case 1 % speed
                    assignin('base', 'wVal', BkptsStruct.w(1));
                    assignin('base', 'frVal', BkptsStruct.fr1(frBkpt));
                 case 2 % fr
                    assignin('base', 'wVal', BkptsStruct.w(wBkpt));
                    assignin('base', 'frVal', BkptsStruct.fr1(1));
             end
             
             out = sim(model);
             
             TnodesExpected = squeeze(out.yout{1}.Values.Data);
             TnodesGridded = squeeze(out.yout{2}.Values.Data)';
             TnodesLPV = squeeze(out.yout{3}.Values.Data)';

             test.verifyEqual(TnodesGridded, TnodesExpected, 'AbsTol', 1e-4, 'RelTol', 1e-4, ...
                 'CASE 1: griddedInterpolant version of InterpolatedStateSpaceThermalModel does not match Simscape baseline');
             test.verifyEqual(TnodesLPV, TnodesExpected, 'AbsTol', 1e-4, 'RelTol', 1e-4, ...
                 'CASE 1: LPV version of InterpolatedStateSpaceThermalModel does not match Simscape baseline');

             % CASE 2: 2nd bkpt -----
             
             assignin('base', 'CapMat', CapMat2);
             assignin('base', 'ResMat', ResMat2);
             switch dimInterp
                 case 1 % speed
                    assignin('base', 'wVal', BkptsStruct.w(2));
                    assignin('base', 'frVal', BkptsStruct.fr1(frBkpt));
                 case 2 % fr
                    assignin('base', 'wVal', BkptsStruct.w(wBkpt));
                    assignin('base', 'frVal', BkptsStruct.fr1(2));
             end

             out = sim(model);
             
             TnodesExpected = squeeze(out.yout{1}.Values.Data);
             TnodesGridded = squeeze(out.yout{2}.Values.Data)';
             TnodesLPV = squeeze(out.yout{3}.Values.Data)';

             test.verifyEqual(TnodesGridded, TnodesExpected, 'AbsTol', 1e-4, 'RelTol', 1e-4, ...
                 'CASE 2: griddedInterpolant version of InterpolatedStateSpaceThermalModel does not match Simscape baseline');
             test.verifyEqual(TnodesLPV, TnodesExpected, 'AbsTol', 1e-4, 'RelTol', 1e-4, ...
                 'CASE 2: LPV version of InterpolatedStateSpaceThermalModel does not match Simscape baseline');

         end
     
     end
end

function samplingGridStruct = getSamplingGridFromBkpts(BkptsStruct)
    % Returns the sampling grid structure based on the BkptsStruct

    fieldNamesCell = fieldnames(BkptsStruct);
    fieldValuesCell = cell(size(fieldNamesCell));
    numFields = length(fieldNamesCell);
    for idxField = 1:numFields
        fieldValuesCell{idxField} = BkptsStruct.(fieldNamesCell{idxField});
    end

    samplingGridValues = cell(size(fieldNamesCell));
    [samplingGridValues{:}] = ndgrid(fieldValuesCell{:});

    samplingGridStruct = cell2struct(samplingGridValues, fieldNamesCell, 1);

end


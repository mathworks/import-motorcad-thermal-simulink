function [Amat, Bmat, ResMat] = thermalStateSpaceInterpolatorImplementation(rotorspeed, flowrate, TcoolantIn, ...
                             Xstruct, AMatND, BMatND, ResMatND, numStates)
    % THERMALSTATESPACEINTERPOLATORIMPLEMENTATION: Interpolates
    % N-dimensional matrices AMatND, BmatND, ResMatND into 2-dimensional
    % matrices Amat, Bmat, ResMat, based on rotorspeed, flowrate, and 
    % TcoolantIn values.
    
    % Copyright 2022 The MathWorks, Inc.
    
    coder.extrinsic('griddedInterpolant')
    
    Amat = zeros(numStates,numStates); 
    Bmat = zeros(numStates,numStates); 
    ResMat = zeros(numStates,numStates);

    Xcell = struct2cell(Xstruct);
    AmatInterpolant = griddedInterpolant(Xcell, AMatND, 'linear', 'linear');
    BmatInterpolant = griddedInterpolant(Xcell, BMatND, 'linear', 'linear');
    ResMatInterpolant = griddedInterpolant(Xcell, ResMatND, 'linear', 'linear');

    Xqcell = num2cell([rotorspeed; flowrate(:); TcoolantIn(:)]);
    Amat(:) = AmatInterpolant(Xqcell);
    Bmat(:) = BmatInterpolant(Xqcell);
    ResMat(:) = ResMatInterpolant(Xqcell);

end




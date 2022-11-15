function [ResMat, CapMat] = getRCFromAB(Amat, Bmat)
    % getRCFromAB Returns the resistance and capacitance matrices corresponding
    % to the A,B state-space matrices.
    
    % Copyright 2022 The MathWorks, Inc.
    
    % Size pre-allocation
    [N,~]=size(Amat);
    ResMat = zeros(N,N);
    C_inv = diag(Bmat);
    CapMat = 1./C_inv;
    
    % Top right half of ResMat matrix
    for row=1:1:N
        for col=row+1:1:N
            ResMat(row,col) = C_inv(row,1)/Amat(row,col);
        end
    end
    
    % Bottom left half of ResMat matrix
    ResMat = ResMat + triu(ResMat,1)';

end

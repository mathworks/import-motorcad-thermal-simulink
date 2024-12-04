function nodelossCorrection = coolantHeatflowAdjustmentImplementation(Tnodes, ResMat, CoolAdjMat, CoolantArrayIdxs, numStates)
    % COOLANTHEATFLOWADJUSTMENTIMPLEMENTATION Calculates node loss correction in upstream coolant nodes.
    % Returns the node loss correction for the coolant nodes that have
    % upstream heat flow. Heat flow can only be conduced downstream, hence
    % any upstream heat flow must be substracted. This function computes
    % this adjustment.

    % Copyright 2022 The MathWorks, Inc.

    nodelossCorrection = zeros(numStates,1);

    for idx1 = 1:length(CoolantArrayIdxs)
        thisCoolNodeIdx = CoolantArrayIdxs(idx1);
        nextCoolNodeIdxs = find(CoolAdjMat(thisCoolNodeIdx,:)==1);
        for idx2 = 1:length(nextCoolNodeIdxs) % for each downstream node
            thisNextCoolNodeIdx = nextCoolNodeIdxs(idx2);       
            ResVal = ResMat(thisCoolNodeIdx, thisNextCoolNodeIdx); % resistance between this upstream and downstream coolant nodes
            powCorrection = max((Tnodes(thisNextCoolNodeIdx) - Tnodes(thisCoolNodeIdx))/ResVal, 0);   % compute heat flow adjustment using node temperature drop and resistance.    
            nodelossCorrection(thisCoolNodeIdx) = nodelossCorrection(thisCoolNodeIdx) - powCorrection;
        end
    end

end




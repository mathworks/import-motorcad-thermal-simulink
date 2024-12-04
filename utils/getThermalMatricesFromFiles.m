function [CapMat, ResMat, PowMat, TempMat, NodeNames] = getThermalMatricesFromFiles(motFileName)
    %GETTHERMALMATRICESFROMFILES Get the capacitance, resistance, and power
    %arrays for thermal model .cmf, .rmf, .pmf text files in the current
    %directory.

    % Copyright 2022-2023 The MathWorks, Inc.

    CapFileName = strrep(motFileName, '.mot', '.cmf');
    ResFileName = strrep(motFileName, '.mot', '.rmf');
    PowFileName = strrep(motFileName, '.mot', '.pmf');
    TempFileName = strrep(motFileName, '.mot', '.tmf');


    % Read cmf file
    [CapMatNm1, ~] = readMfFile1D(CapFileName); % size numNodes-1
    numNodes = length(CapMatNm1) + 1; % number of nodes (+1 for ambient)
    CapMat = zeros(numNodes,1);
    CapMat(2:end) = CapMatNm1;

    % Read tmf file
    [TempMat, NodeNames] = readMfFile1D(TempFileName); % size numNodes
    
    CapMat(1,1) = 1e20; % Ambient must have large capacitance
    for idx_node=2:numNodes      
        if  TempMat(idx_node,1) ~= -10000000 % Fixed-temperature node
            % Model these nodes as having very large capcitances
            CapMat(idx_node,1) = 1e20;
        end 
    end

    % Read rmf file
    ResMat = readMfFile2D(ResFileName, numNodes);
    
    % Read pmf file
    [PowMatNm1,~] = readMfFile1D(PowFileName); % size numNodes-1
    PowMat = zeros(numNodes,1);
    PowMat(2:end) = PowMatNm1;

end

function [array1d, nodeNames] = readMfFile1D(FileName)

    fid = fopen(FileName, 'rt');
    fgetl(fid); % first line not useful
    fgetl(fid); % second line contains number of nodes
    
    datacell = textscan(fid, '%d%s%f', 'Delimiter',';', 'CollectOutput', 1);

    nodeNames = datacell{2};     %a cell array of strings
    array1d = datacell{3};    %as a numeric array

    fclose(fid);
end

function array2d = readMfFile2D(FileName, numNodes)

    fileID = fopen(FileName,'r');

    % If error exit
    if fileID < 0
        disp('Error loading the file, exit')
        return
    end
    
    % Read first lines without useful info
    textscan(fileID, '%[^\n\r]', 3, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    % Read resistance values
    formatSpec = '%s';    
    % Read the resistance values. We need one more node for the ambient now!!!
    for i=1:numNodes
        formatSpec = [formatSpec, '%f'];
    end
    formatSpec = [formatSpec, '%[^\n\r]'];
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ';', 'TextType', 'string', 'ReturnOnError', false);
    fclose(fileID);
    
    % Generate the resistace matrix
    array2d = zeros(numNodes);
    
    % Create the matrix (R matrix in only 1 vector)
    dtmp = cell2mat(dataArray(1,2:numNodes+1));
    array2d(:) = dtmp(:);
end
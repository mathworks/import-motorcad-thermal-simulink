function [GroupNamesAndMcadIdxes, NodeNamesAndMcadIdx] = getGroupNamesAndNodes(motFileName)
    %GETGROUPNAMESANDNODES Get the group names, node names, and node
    %MotorCAD index numbers from a .mot file. Requires thermal matrices
    %previously exported.
    % Input arguments:
    % - motFileName: [string/char]: file name of a .mot file (including .mot)
    % Output arguments:
    % - GroupNames: [Nx2 cell]: cell array with the groups information.
    % Each row represents a group. The first column is the group name, the
    % second column is an array with the Mcad indices of the nodes forming
    % this group.
    % - NodeNamesAndMcadIdx: [Mx2 cell]: cell array with the nodes
    % information. Each row represents a node. The first column is the node
    % name, the second column is the Mcad index associated to the node.

    % Copyright 2022 The MathWorks, Inc.

    NodeFileName = strrep(motFileName, '.mot', '.nmf');

    % Read nmf file
    [GroupNames, NodeNames, GroupIdxs, NodeIdx] = readNmfFile(NodeFileName); % size numNodes-1

    % remove initial and final parenthesis from obj.NodeNames
    for idxNode = 1:length(NodeNames)
        NodeNames{idxNode} = NodeNames{idxNode}(2:end-1);
    end
    

    numNodes = length(NodeIdx);
    numGroups = length(GroupNames);

    GroupNamesAndMcadIdxes = cell(numGroups, 2);
    GroupNamesAndMcadIdxes(:,1) = GroupNames;
    GroupNamesAndMcadIdxes(:,2) = GroupIdxs;

    NodeNamesAndMcadIdx = cell(numNodes, 2);
    NodeNamesAndMcadIdx(:,1) = NodeNames;
    NodeNamesAndMcadIdx(:,2) = NodeIdx;
   
end

function [GroupNames, NodeNames, GroupIdxs, NodeIdx] = readNmfFile(FileName)

    fid = fopen(FileName, 'rt');
    fgetl(fid); % first line not useful
    fgetl(fid); % second line contains number of nodes
    fgetl(fid); % third line is whitespace

    datacell = textscan(fid, '%s', 'CollectOutput', 1, 'Delimiter', newline);
    datacell = datacell{1}; % data is one level below

    numGroups = 0;
    numNodes = 0;
    GroupNames = {};
    NodeNames = {};
    GroupIdxs = {};
    NodeIdx = {};
    for idxLine = 1:length(datacell)
        thisLine = datacell{idxLine};
        if ~isempty(thisLine) % skip empty lines
            if startsWith(thisLine, '[') % new group
                % New group           
                numGroups = numGroups+1;
                GroupNames{end+1} = thisLine(2:end-1);     %#ok<AGROW> 
                GroupIdxs{end+1} = []; %#ok<AGROW> 
            else
                % Keep adding to the current group
                numNodes = numNodes+1;
                nodeIdx = str2double(extractBefore(thisLine, ' '));
                nodeName = extractAfter(thisLine, '('); % Assumes node name between parenthesis
                nodeName = nodeName(1:end-1); % last character is ")"
                NodeNames{end+1} = nodeName; %#ok<AGROW> 
                NodeIdx{end+1} = nodeIdx; %#ok<AGROW> 
                GroupIdxs{end}(end+1) = nodeIdx;
            end
        end
    end

    fclose(fid);
end

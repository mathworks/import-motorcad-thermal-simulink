classdef ThermalInterface < mcadinterface.BasicInterface
    %THERMALINTERFACE Motor-CAD interface for thermal modeling workflows.
   
    % Copyright 2022-2025 The MathWorks, Inc.

    properties(SetAccess=protected)
        workingDirectory % Working directory
        mcadROMLibName % Library of Simulink subsystems required for the ROM model
    end

    properties(SetAccess=private)
        NodeNames (:,1) cell % List of node names
        NodeNamesAndMcadIdx (:,2) cell % List of node names and the corresponding Motor-CAD index
        GroupNamesAndMcadIdxes (:,2) cell % List of group names and the Motor-CAD index of the nodes in each group
        CoolingSystemNamesAndMcadIdxes (:,2) cell % List of cooling system group names and the Motor-CAD index of the fluid nodes in each cooling system

        % Lab Maps

        Speed_Mat (:,:) double % Speed map [rpm]
        Shaft_Torque_Mat (:,:) double % Torque map [N*m]
        Stator_Copper_Loss_Mat (:,:) double % Stator copper loss map [W]
        Rotor_Cage_Loss_Mat (:,:) double % Rotor copper loss map [W]
        Iron_Loss_Stator_Back_Iron_Mat (:,:) double % Stator back iron loss map [W]
        Iron_Loss_Stator_Tooth_Mat (:,:) double % Stator tooth iron loss map [W]
        Stray_Load_Loss_Mat (:,:) double % Stray loss map [W]
        Magnet_Loss_Mat (:,:) double % Magnet loss map [W]
        Iron_Loss_Rotor_Pole_Mat (:,:) double % Rotor pole iron loss map [W]
        Iron_Loss_Rotor_Back_Iron_Mat (:,:) double % Rotor back iron loss map [W]
        Iron_Loss_Rotor_Tooth_Mat (:,:) double % Rotor tooth iron loss map [W]
        Friction_Loss_Mat (:,:) double % Friction loss map [W]
        Windage_Loss_Mat (:,:) double % Windage loss map [W]
        Stator_Copper_Loss_AC_Mat (:,:) double % Stator copper AC loss map [W]
        Banding_Loss_Mat (:,:) double % Banding loss map [W]
        Sleeve_Loss_Mat (:,:) double % Windage loss map [W]

        % Thermal matrices

        CapMat (:,:) double % Node capacitance vector
        ResMat (:,:) double % Node resistance matrix
        PowMat (:,:) double % Node steady-state power vector
        TempMat (:,:) double % Node temperature boundary condition

        % State-space matrices

        Amat (:,:) double % State-space A matrix
        Bmat (:,:) double % State-space B matrix

        % Cooling system data      

        EnabledCoolingSystems (:,1) logical % List of which cooling systems are enabled (true or false)
        CoolingSystemsDigraphs (:,1) struct % Structure containing the directed graph of each cooling systems node flow
        AdjacencyMat (:,:) double % Global adjacency matrix for cooling systems node flow connectivity
        InletArrayIdxs (:,1) double % Inlet nodes indices
        OutletArrayIdxs (:,1) double % Outlet nodes indices
        CoolantArrayIdxs (:,1) double % Coolant nodes indices
    end

    properties(Constant)
        % List of supported cooling systems
        SupportedCoolingSystemsMcadNames = ...
                    {'Ventilated'; ...
                     'Housing Water Jacket'; ...
                     'Shaft Spiral Groove'; ...
                     'Wet Rotor'; ...
                     'Spray Cooling'; ...
                     'Rotor Water Jacket'; ...
                     'Slot Water Jacket'; ...
                                   }; 
    end

    methods(Access=public)

        function obj = ThermalInterface(motFile)
            % ThermalInterface constructor
            obj = obj@mcadinterface.BasicInterface(motFile);
            fullFilePath = which(motFile);
            obj.workingDirectory = fileparts(fullFilePath); % use same folder as .mot file
            obj.mcadROMLibName = 'mcadROM_lib';

            obj.calculateThermalSteadyState(); % initialize node temperatures
            obj.updateModel();
            
        end

        function updateModel(obj)
            % Update state-space matrices, group names, nodes, cooling
            % systems, and loss tables, based on the current state of
            % Motor-CAD variables.

            % Set thermal matrices, group names, node names, and mcad idxs
            obj.updateMatricesAndGroupNamesAndNodes();

            obj.updateCoolingSystemsData();
            obj.checkCoolingSystemsIndependent();

            obj.updateLabLossTables();

        end

        % Basic MotorCAD thermal commands ---------
        function calculateLabOperatingPoint(obj)
            % Calculate Lab Operating Point
            obj.mcad.calculate_operating_point_lab();
        end

        function calculateThermalSteadyState(obj)
            % Calculate thermal steady state
            obj.mcad.do_steady_state_analysis();
        end

        function calculateThermalTransient(obj)
            % Calculate thermal transient
            obj.mcad.do_transient_analysis();
        end

        function runThermalSteadyStateWithSpecifiedLosses(obj, lossVec)
            % Run a thermal steady-state calculation with specified loss
            % values for each type of loss in "LossTypes".
            
            obj.ThermalSteadyOrTransientChoice = 0;
            obj.MagneticThermalCouplingChoice = 0;
            obj.LabThermalCouplingChoice = 0;
            obj.LossValues = lossVec;

            obj.calculateThermalSteadyState();

        end

        function runThermalSteadyStateWithSpecifiedTorqueSpeed(obj, torqueVal, speedVal)
            % Run a thermal steady-state calculation for a specified operating point (load
            % torque and shaft speed).

            obj.ThermalSteadyOrTransientChoice = 0;
            obj.MagneticThermalCouplingChoice = 0;
            obj.LabThermalCouplingChoice = 2;
            obj.Shaft_Torque_Nm = torqueVal;
            obj.Shaft_Speed_RPM = speedVal;

            obj.calculateLabOperatingPoint();

        end

        function runThermalTransientWithSpecifiedLosses(obj, lossVec, stopTime, numTimeSteps)
            % Run a thermal transient calculation with specified loss
            % values for each type of loss in "LossTypes". You also specify the stop
            % time and number of time steps

            obj.ThermalSteadyOrTransientChoice = 1;
            obj.TransientOption = 0;           
            obj.StopTimeChoice = 0; % fixed period
            obj.StopTime = stopTime;
            obj.TransientLossesOrTorqueChoice = 0; % specify losses
            obj.NumTimeSteps = numTimeSteps;
            obj.InitialTemperatureOption = 0; % initialize at ambient temperature

            obj.LossValues = lossVec;

            obj.calculateThermalTransient();

        end

        function runThermalTransientWithSpecifiedTorqueSpeed(obj, torqueVal, speedVal, stopTime, numTimeSteps)
            % Run a thermal transient calculation with specified operating point (load
            % torque and shaft speed). You also specify the stop time and
            % number of time steps.

            obj.ThermalSteadyOrTransientChoice = 1;
            obj.TransientOption = 0;           
            obj.StopTimeChoice = 0; % fixed period
            obj.StopTime = stopTime;
            obj.TransientLossesOrTorqueChoice = 1; % specify torque
            obj.NumTimeSteps = numTimeSteps;
            obj.InitialTemperatureOption = 0; % initialize at ambient temperature

            obj.TransientTorqueValue = torqueVal;
            obj.Shaft_Speed_RPM = speedVal;

            obj.calculateThermalTransient();

        end

        function writeThermalStateSpaceFiles(obj)
            % Export thermal matrices into text files (.cmf, .rmf, .pmf,
            % .tmf, .nmf)
            obj.mcad.export_matrices(obj.workingDirectory);
        end

        function [CapMat, ResMat, PowMat, TempMat, NodeNames] = getThermalStateSpaceMatricesFromFiles(obj)
            %Get the capacitance, resistance, and power
            %arrays for thermal model .cmf, .rmf, .pmf text files in the current
            %directory.

            motFileName = obj.motFullFile;

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
            % remove initial and final parenthesis from obj.NodeNames
            for idxNode = 1:length(NodeNames)
                NodeNames{idxNode} = NodeNames{idxNode}(2:end-1);
            end
            
            CapMat(1,1) = 1e20; % Ambient must have large capacitance
            for idx_node=2:numNodes      
                if  TempMat(idx_node,1) ~= -10000000 % Fixed-temperature node
                    % Model these nodes as having very large capacitances
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

        function [GroupNamesAndMcadIdxes, NodeNamesAndMcadIdx] = updateMatricesAndGroupNamesAndNodes(obj)
            % Update state-space matrices, group data (names and Motor-CAD
            % indices) and node data (names and Motor-CAD indices)

            xxLossValues = obj.LossValues;
            xxEnableStatorTempCoeffRes = obj.EnableStatorTempCoeffRes;
            xxEnableRotorTempCoeffRes = obj.EnableRotorTempCoeffRes;

            obj.setupForThermalMatricesUpdate();
            obj.writeThermalStateSpaceFiles();
            [xCapMat, xResMat, xPowMat, xTempMat, xNodeNames] = obj.getThermalStateSpaceMatricesFromFiles();

            % restore original properties prior to matrices update
            obj.LossValues = xxLossValues; 
            obj.EnableStatorTempCoeffRes = xxEnableStatorTempCoeffRes;
            obj.EnableRotorTempCoeffRes = xxEnableRotorTempCoeffRes;

            obj.CapMat = xCapMat;
            obj.ResMat = xResMat;
            obj.PowMat = xPowMat;
            obj.TempMat = xTempMat;  

            [xAmat, xBmat] = getStateSpaceMatricesFromThermalMatrices(xCapMat, xResMat);
            obj.Amat = xAmat;
            obj.Bmat = xBmat;

            obj.NodeNames = xNodeNames;
        
            % Read nmf file
            motFileName = obj.motFullFile;
            NodeFileName = strrep(motFileName, '.mot', '.nmf');        
            
            [GroupNames, nodeNames, GroupIdxs, NodeIdx] = readNmfFile(NodeFileName); % size numNodes-1      
        
            numNodes = length(NodeIdx);
            numGroups = length(GroupNames);
        
            GroupNamesAndMcadIdxes = cell(numGroups, 2);
            GroupNamesAndMcadIdxes(:,1) = GroupNames;
            GroupNamesAndMcadIdxes(:,2) = GroupIdxs;
        
            NodeNamesAndMcadIdx = cell(numNodes, 2);
            NodeNamesAndMcadIdx(:,1) = nodeNames;  
            NodeNamesAndMcadIdx(:,2) = NodeIdx;

            % sort according to obj.NodeNames
            sortIdxs = nan(size(nodeNames));
            for idxNode = 1:numNodes                
                sortIdxs(idxNode) = find(strcmp(NodeNamesAndMcadIdx(:,1),obj.NodeNames{idxNode}));
            end
            NodeNamesAndMcadIdx(:,2) = NodeNamesAndMcadIdx(sortIdxs,2);        
            for idxNode = 1:numNodes   
                NodeNamesAndMcadIdx(idxNode,1) = obj.NodeNames(idxNode);
            end

            % remove initial and final parenthesis from obj.NodeNames
            obj.NodeNames = NodeNamesAndMcadIdx(:,1);

            % Update
            obj.GroupNamesAndMcadIdxes = GroupNamesAndMcadIdxes;
            obj.NodeNamesAndMcadIdx = NodeNamesAndMcadIdx;

        end

        function updateCoolingSystemsData(obj)
            % Update cooling systems data (AdjacencyMat, InletArrayIdxs,
            % OutletArrayIdxs, CoolantArrayIdxs)

            groupNames = obj.GroupNamesAndMcadIdxes(:,1);
            obj.CoolingSystemNamesAndMcadIdxes = obj.GroupNamesAndMcadIdxes(contains(groupNames, obj.SupportedCoolingSystemsMcadNames), :);
            obj.EnabledCoolingSystems = contains(obj.SupportedCoolingSystemsMcadNames, obj.CoolingSystemNamesAndMcadIdxes(:,1));
            numSuppCoolSys = length(obj.SupportedCoolingSystemsMcadNames);
            obj.CoolingSystemsDigraphs = struct();
            for idxSuppCoolSys = 1:numSuppCoolSys
                isThisCoolSysEnabled = obj.EnabledCoolingSystems(idxSuppCoolSys);
                thisCoolantGroupName = obj.SupportedCoolingSystemsMcadNames{idxSuppCoolSys};
                if isThisCoolSysEnabled
                    [thisCoolantDigraph, thisCoolantNodeNamesAndArrayIdx] = obj.getDigraphForCoolantGroup(thisCoolantGroupName, false); % don't draw plot
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).Digraph = thisCoolantDigraph;
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).NodeNamesAndArrayIdx = thisCoolantNodeNamesAndArrayIdx;
                else
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).Digraph = []; % empty digraph
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).NodeNamesAndArrayIdx = []; % empty NodeNamesAndArrayIdx
                end
            end

            [xAdjacencyMat, xInletArrayIdxs, xOutletArrayIdxs, xCoolantArrayIdxs] = getAdjacencyMatAndInletOutletIdxs(obj);
            obj.AdjacencyMat = xAdjacencyMat;
            obj.InletArrayIdxs = xInletArrayIdxs;
            obj.OutletArrayIdxs = xOutletArrayIdxs;
            obj.CoolantArrayIdxs = xCoolantArrayIdxs;

        end

        function TnodesVec = getSteadyStateTemperatureForNodeMcadIdxs(obj, mcadIdxs) 
            % Get steady-state temperature for a list of nodes specified by
            % the Motor-CAD indices.

            TnodesVec = nan(length(mcadIdxs), 1);
            for idxNode = 1:length(mcadIdxs)
                Tval = obj.mcad.get_node_temperature(int32(mcadIdxs(idxNode)));
                TnodesVec(idxNode) = double(Tval);
            end
        end

        function PnodesVec = getSteadyStatePowerForNodeMcadIdxs(obj, mcadIdxs) 
            % Get steady-state power for a list of nodes specified by
            % the Motor-CAD indices.

            PnodesVec = nan(length(mcadIdxs), 1);
            for idxNode = 1:length(mcadIdxs)
                Pval = obj.mcad.get_node_power(int32(mcadIdxs(idxNode)));
                PnodesVec(idxNode) = double(Pval);

            end
        end

        function [tVec, TnodesVec] = getTransientTemperatureForNodeMcadIdxs(obj, mcadIdxs)
            % Get transient temperature time series for a list of nodes specified by
            % the Motor-CAD indices.

            numTimePoints = double(obj.mcad.get_variable('Simple_Transient_Number_Points'));
            numTimePoints = numTimePoints+1; % include initial time
            TnodesVec = nan(length(mcadIdxs), numTimePoints);
            tVec = nan(1,numTimePoints);
            for idxNode = 1:length(mcadIdxs)
                for idxTime = 1:numTimePoints
                    tuple = obj.mcad.get_temperature_graph_point(int32(mcadIdxs(idxNode)),int32(idxTime-1));
                    tuple = cell(tuple);
                    x = double(tuple{1});
                    y = double(tuple{2});
                    tVec(idxTime) = x;
                    TnodesVec(idxNode, idxTime) = y;
                end
            end
        end

        function [tVec, PnodesVec] = getTransientPowerForNodeMcadIdxs(obj, mcadIdxs)
            % Get transient power time series for a list of nodes specified by
            % the Motor-CAD indices.

            numTimePoints = double(obj.mcad.get_variable('Simple_Transient_Number_Points'));
            numTimePoints = numTimePoints+1; % include initial time
            PnodesVec = nan(length(mcadIdxs), numTimePoints);
            tVec = nan(1,numTimePoints);
            for idxNode = 1:length(mcadIdxs)
                for idxTime = 1:numTimePoints
                    tuple = obj.mcad.get_power_graph_point(int32(mcadIdxs(idxNode)),int32(idxTime-1));
                    tuple = cell(tuple);
                    x = double(tuple{1});
                    y = double(tuple{2});
                    tVec(idxTime) = x;
                    PnodesVec(idxNode, idxTime) = y;
                end
            end
            PnodesVec(1,:) = 0; 
        end

        function plotTransientTemperatureForNodeMcadIdxs(obj, mcadIdxs)
            % Plot transient temperature time series for a list of nodes specified by
            % the Motor-CAD indices.

            [tVec, TnodesVec] = obj.getTransientTemperatureForNodeMcadIdxs(mcadIdxs);
            allMcadIdxs = [obj.NodeNamesAndMcadIdx{:,2}];
            legendNodeNames = cell(length(mcadIdxs),1);
            for idxNode = 1:length(mcadIdxs)
                legendNodeNames{idxNode} = obj.NodeNames{allMcadIdxs == mcadIdxs(idxNode)};
            end
            figure();
            plot(tVec, TnodesVec);
            legend(legendNodeNames, 'Interpreter', 'none');
            grid on
            xlabel('Time [s]');
            ylabel('Node temperature [degC]');

        end

        % Basic MotorCAD Lab commands --------------
        function calculateMagneticLab(obj)
            % Calculate magnetic lab

            obj.mcad.calculate_magnetic_lab()
            % calculate_magnetic_lab uses the existing pre-built magnetic
            % lab (from a previous call to BuildModel_Lab)
        end

        function updateLabLossTables(obj)
            % Update Lab loss tables for each type of loss

            [~, motFile] = fileparts(obj.motFullFile);
            LabFileName = fullfile(extractBefore(obj.motFullFile, motFile), motFile, 'Lab', 'MotorLAB_elecdata.mat');
            if ~isfile(LabFileName) % need to generate Lab mat file
                obj.calculateMagneticLab(); 
            end
            outLoad = load(LabFileName);
        
            % Speed and Shaft_Torque always exist
            Speed = outLoad.Speed;
            Shaft_Torque = outLoad.Shaft_Torque;
            % Losses may be specific to each machine - check first
            if isfield(outLoad, 'Stator_Copper_Loss')
                Stator_Copper_Loss = outLoad.Stator_Copper_Loss;
            else % this machine type does not include this loss type
                Stator_Copper_Loss = zeros(size(Speed));
            end
            if isfield(outLoad, 'Rotor_Cage_Loss')
                Rotor_Cage_Loss = outLoad.Rotor_Cage_Loss;
            else % this machine type does not include this loss type
                Rotor_Cage_Loss = zeros(size(Speed));
            end
            if isfield(outLoad, 'Iron_Loss_Stator_Back_Iron')
                Iron_Loss_Stator_Back_Iron = outLoad.Iron_Loss_Stator_Back_Iron;
            else % this machine type does not include this loss type
                Iron_Loss_Stator_Back_Iron = zeros(size(Speed));
            end
            if isfield(outLoad, 'Iron_Loss_Stator_Tooth')
                Iron_Loss_Stator_Tooth = outLoad.Iron_Loss_Stator_Tooth;
            else % this machine type does not include this loss type
                Iron_Loss_Stator_Tooth = zeros(size(Speed));
            end
            if isfield(outLoad, 'Stray_Load_Loss')
                Stray_Load_Loss = outLoad.Stray_Load_Loss;
            else % this machine type does not include this loss type
                Stray_Load_Loss = zeros(size(Speed));
            end
            if isfield(outLoad, 'Magnet_Loss')
                Magnet_Loss = outLoad.Magnet_Loss;
            else % this machine type does not include this loss type
                Magnet_Loss = zeros(size(Speed));
            end
            if isfield(outLoad, 'Iron_Loss_Rotor_Pole')
                Iron_Loss_Rotor_Pole = outLoad.Iron_Loss_Rotor_Pole;
            else % this machine type does not include this loss type
                Iron_Loss_Rotor_Pole = zeros(size(Speed));
            end
            if isfield(outLoad, 'Iron_Loss_Rotor_Back_Iron')
                Iron_Loss_Rotor_Back_Iron = outLoad.Iron_Loss_Rotor_Back_Iron;
            else % this machine type does not include this loss type
                Iron_Loss_Rotor_Back_Iron = zeros(size(Speed));
            end
            if isfield(outLoad, 'Iron_Loss_Rotor_Tooth')
                Iron_Loss_Rotor_Tooth = outLoad.Iron_Loss_Rotor_Tooth;
            else % this machine type does not include this loss type
                Iron_Loss_Rotor_Tooth = zeros(size(Speed));
            end
            if isfield(outLoad, 'Friction_Loss')
                Friction_Loss = outLoad.Friction_Loss;
            else % this machine type does not include this loss type
                Friction_Loss = zeros(size(Speed));
            end
            if isfield(outLoad, 'Windage_Loss')
                Windage_Loss = outLoad.Windage_Loss;
            else % this machine type does not include this loss type
                Windage_Loss = zeros(size(Speed));
            end
            if isfield(outLoad, 'Stator_Copper_Loss_AC')
                Stator_Copper_Loss_AC = outLoad.Stator_Copper_Loss_AC;
            else % this machine type does not include this loss type
                Stator_Copper_Loss_AC = zeros(size(Speed));
            end
            if isfield(outLoad, 'Banding_Loss')
                Banding_Loss = outLoad.Banding_Loss;
            else % this machine type does not include this loss type
                Banding_Loss = zeros(size(Speed));
            end
            if isfield(outLoad, 'Sleeve_Loss')
                Sleeve_Loss = outLoad.Sleeve_Loss;
            else % this machine type does not include this loss type
                Sleeve_Loss = zeros(size(Speed));
            end

            obj.Speed_Mat = Speed;
            obj.Shaft_Torque_Mat = Shaft_Torque;
            obj.Stator_Copper_Loss_Mat = Stator_Copper_Loss;
            obj.Rotor_Cage_Loss_Mat = Rotor_Cage_Loss;
            obj.Iron_Loss_Stator_Back_Iron_Mat = Iron_Loss_Stator_Back_Iron;
            obj.Iron_Loss_Stator_Tooth_Mat = Iron_Loss_Stator_Tooth;
            obj.Stray_Load_Loss_Mat = Stray_Load_Loss;
            obj.Magnet_Loss_Mat = Magnet_Loss;
            obj.Iron_Loss_Rotor_Pole_Mat = Iron_Loss_Rotor_Pole;
            obj.Iron_Loss_Rotor_Back_Iron_Mat = Iron_Loss_Rotor_Back_Iron;
            obj.Iron_Loss_Rotor_Tooth_Mat = Iron_Loss_Rotor_Tooth;
            obj.Friction_Loss_Mat = Friction_Loss;
            obj.Windage_Loss_Mat = Windage_Loss;
            obj.Stator_Copper_Loss_AC_Mat = Stator_Copper_Loss_AC;
            obj.Banding_Loss_Mat = Banding_Loss;
            obj.Sleeve_Loss_Mat = Sleeve_Loss;
    
        end

        % Advanced workflows -----------

        function [DirectedGraphCoolant, NodeNamesAndArrayIdx] = getDigraphForCoolantGroup(obj, CoolantGroupName, plotCoolantGraphFlag)
            % Get the directed flow graph for a specific cooling system
            % group. Optionally plot the graph.
            % Input arguments:
            % - CoolantGroupName [string] : name of the cooling group
            % Output arguments:
            % - Digraph [1x1 digraph] : Directed graph that represents the coolant flow from inlet to
            % outlet nodes with possible path splits and merges
            % - NodeNamesAndArrayIdx [Mx3 cell] : cell with the node name,
            % array index, and Mcad index for each of the coolant nodes. The order is
            % in the same order of the AdjacencyMatrix rows

            % Extract coolant variables
            AllNodeNames = obj.NodeNames;
            assert(all(strcmp(AllNodeNames, obj.NodeNamesAndMcadIdx(:,1)))) % required condition
            AllMcadIdxes = [obj.NodeNamesAndMcadIdx{:,2}];
            xResMat = obj.ResMat;  
            CoolantNodesMcadIdx = obj.CoolingSystemNamesAndMcadIdxes{contains(obj.CoolingSystemNamesAndMcadIdxes(:,1), CoolantGroupName), 2};
            CoolantNodeNames = cell(size(CoolantNodesMcadIdx));
            CoolantArrayIdxes = nan(size(CoolantNodesMcadIdx));
            NodeNamesAndArrayIdx = cell(length(CoolantNodesMcadIdx),3);
            for idxCoolantNode = 1:length(CoolantNodesMcadIdx)
                CoolantNodeNames{idxCoolantNode} = AllNodeNames{CoolantNodesMcadIdx(idxCoolantNode) == AllMcadIdxes};
                CoolantArrayIdxes(idxCoolantNode) = find(CoolantNodesMcadIdx(idxCoolantNode) == AllMcadIdxes);
                NodeNamesAndArrayIdx{idxCoolantNode, 1} = CoolantNodeNames{idxCoolantNode};
                NodeNamesAndArrayIdx{idxCoolantNode, 2} = CoolantArrayIdxes(idxCoolantNode);
                NodeNamesAndArrayIdx{idxCoolantNode, 3} = CoolantNodesMcadIdx(idxCoolantNode);
            end

            ResMatCoolant = xResMat(CoolantArrayIdxes, CoolantArrayIdxes);

            % Inlet to outlet graph search
            isInlet = contains(CoolantNodeNames, 'Inlet');
            if any(isInlet)
                if sum(int32(isInlet)) > 1
                    error('Only one inlet per cooling system is supported');
                else
                    % BFS GRAPH ALGORITHM
                    UndirectedAdjacencyMatrix = abs(ResMatCoolant)<1e8 & ResMatCoolant~=0; % nodes are connected and is not itself
                    GraphCoolant = graph(UndirectedAdjacencyMatrix, CoolantNodeNames);
                    nodeVisitList = GraphCoolant.bfsearch(CoolantNodeNames(isInlet));
                    DirectedGraphCoolant = digraph(adjacency(GraphCoolant), CoolantNodeNames); % all connections are double-edged
                    DirectedGraphCoolant = DirectedGraphCoolant.rmedge(nodeVisitList(2:end), nodeVisitList(1:end-1)); % remove edges from BFS
                    % remove cycles
                    [~,edge_indices] = DirectedGraphCoolant.dfsearch(CoolantNodeNames(isInlet), 'edgetodiscovered', 'Restart', true);
                    DirectedGraphCoolant = DirectedGraphCoolant.rmedge(edge_indices);
                    
                    if plotCoolantGraphFlag
                        figure('Name', CoolantGroupName);
                        graphPlot = DirectedGraphCoolant.plot();
                        graphPlot.Interpreter = 'none'; % avoid subscript labels
                        set(graphPlot.Parent.Title, 'String', CoolantGroupName)
                    end

                end
            else
               error('No inlet found. Cooling groups must contain an inlet node, identified with an "Inlet" substring in the node name.');
            end

        end

        function areIndependent = checkCoolingSystemsIndependent(obj)
            % Returns true if the enabled cooling systems are independent
            % i.e. they don't merge or split into each other (no fluid
            % is mixed between cooling systems)

            % Get ArrayIdx for both of the cooling system nodes
            enabledCoolSysIdxes = find(obj.EnabledCoolingSystems);
            if length(enabledCoolSysIdxes)<=1 
                areIndependent = true;
                return
            end
            pairCasesMat = nchoosek(enabledCoolSysIdxes, 2);
            [numCases, ~] = size(pairCasesMat);
            couplingMat = zeros(length(obj.EnabledCoolingSystems));
            for idxCase = 1:numCases
                firstCoolSys = pairCasesMat(idxCase, 1);
                secondCoolSys = pairCasesMat(idxCase, 2);
                firstCoolStruct = obj.CoolingSystemsDigraphs(firstCoolSys);
                firstCoolArrayIdxs = [firstCoolStruct.NodeNamesAndArrayIdx{:,2}]';
                secondCoolStruct = obj.CoolingSystemsDigraphs(secondCoolSys);
                secondCoolArrayIdxs = [secondCoolStruct.NodeNamesAndArrayIdx{:,2}]';
                CoolantArrayIdxes = cat(1, firstCoolArrayIdxs, secondCoolArrayIdxs);

                % Get mutual adjacency matrix
                xResMat = obj.ResMat;
                ResMatCoolant = xResMat(CoolantArrayIdxes, CoolantArrayIdxes);
                UndirectedAdjacencyMatrix = abs(ResMatCoolant)<1e8 & ResMatCoolant~=0 ;

                % Get non-diagonal submatrice (cross-coupling)
                NondiagonalSubmatrix = UndirectedAdjacencyMatrix(1:length(firstCoolArrayIdxs), (length(firstCoolArrayIdxs)+1):end);
                if any(NondiagonalSubmatrix(:))
                    couplingMat(firstCoolSys, secondCoolSys) = 1;
                    couplingMat(secondCoolSys, firstCoolSys) = 1;
                end

            end

            if any(couplingMat(:)==1)
                areIndependent = false;
                warnString = 'INVALID COOLING SYSTEMS SETTINGS: At least two cooling systems are coupled. This is not supported yet. Modify the Motor-CAD cooling system options to resolve this issue. ';
                warning(warnString);
                warnString = strcat(warnString, [newline, 'Invalid couplings between cooling systems found: ']);
                for idx1 = 1:length(enabledCoolSysIdxes)
                    for idx2 = 1:length(enabledCoolSysIdxes)
                        if couplingMat(idx1,idx2)==1
                            name1 = obj.CoolingSystemNamesAndMcadIdxes{idx1, 1};
                            name2 = obj.CoolingSystemNamesAndMcadIdxes{idx2, 1};
                            warnString = strcat(warnString, [char(name1), ' coupled with ', char(name2), '. ', newline]);
                        end
                    end
                end
                warndlg(warnString); % pop-up warning dialog
            else
                areIndependent = true;
            end

        end

        function [AdjacencyMat, InletArrayIdxs, OutletArrayIdxs, CoolantArrayIdxs] = getAdjacencyMatAndInletOutletIdxs(obj)
            % Get adjacency matrix, inlet node indices, outlet node indices
            % and all coolant node indices for each cooling system.

            IdxsEnabledCoolSys = find(obj.EnabledCoolingSystems);
            InletArrayIdxs = nan(length(IdxsEnabledCoolSys),1);
            OutletArrayIdxs = nan(length(IdxsEnabledCoolSys),1);
            CoolantArrayIdxs = [];
            AdjacencyMat = zeros(size(obj.ResMat));
            for idx = 1:length(IdxsEnabledCoolSys)
                idxCoolSys = IdxsEnabledCoolSys(idx);
                thisDigraph = obj.CoolingSystemsDigraphs(idxCoolSys).Digraph;
                thisNodeNamesAndArrayIdx = obj.CoolingSystemsDigraphs(idxCoolSys).NodeNamesAndArrayIdx;
                nodeNames = thisNodeNamesAndArrayIdx(:,1);
                nodeArrayIdxs = [thisNodeNamesAndArrayIdx{:,2}];
                CoolantArrayIdxs = cat(1,CoolantArrayIdxs, nodeArrayIdxs(:));
                % find inlet and outlet
                inletIdx = contains(nodeNames, 'Inlet');
                pathSequence = thisDigraph.bfsearch(find(inletIdx));
                outletIdx = pathSequence(end);
                InletArrayIdxs(idx) = nodeArrayIdxs(inletIdx);
                OutletArrayIdxs(idx) = nodeArrayIdxs(outletIdx);
                % Add edges to global adjacency matrix
                AdjacencyMat(nodeArrayIdxs, nodeArrayIdxs) = adjacency(thisDigraph);
            end
                
        end

        function lossDistrForEachType = getLossDistrForEachType(obj)
            % Get loss distribution amongst nodes, for each loss type.

            xxLossValues = obj.LossValues;
            xxEnableStatorTempCoeffRes = obj.EnableStatorTempCoeffRes;
            xxEnableRotorTempCoeffRes = obj.EnableRotorTempCoeffRes;

            obj.EnableStatorTempCoeffRes = 0;
            obj.EnableRotorTempCoeffRes = 0;
            numLossTypes = length(obj.LossValues);
            numStates = length(obj.NodeNames);
            totalLossForTest = 100; % [W]
            lossDistrForEachType = zeros(numLossTypes, numStates);
            for idxLossType = 1:numLossTypes
                lossVals = zeros(numLossTypes,1);
                lossVals(idxLossType) = totalLossForTest;
                obj.runThermalSteadyStateWithSpecifiedLosses(lossVals);
                obj.writeThermalStateSpaceFiles();
                [~, ~, xPowMat, ~, ~] = obj.getThermalStateSpaceMatricesFromFiles();
                xPowMat(xPowMat<0)=0;
                powLossDistribution = xPowMat./totalLossForTest;
                if (sum(powLossDistribution)>1e-3) && (abs(sum(powLossDistribution)-1)>1e-3)
                    warning(['The loss distribution for ' char(obj.LossTypes{idxLossType}), ' is not consistent. Results may be inaccurate.']);
                end
                lossDistrForEachType(idxLossType,:) = powLossDistribution;
            end
            obj.EnableStatorTempCoeffRes = 1;
            obj.EnableRotorTempCoeffRes = 1;

            % restore original properties prior to matrices update
            obj.LossValues = xxLossValues; 
            obj.EnableStatorTempCoeffRes = xxEnableStatorTempCoeffRes;
            obj.EnableRotorTempCoeffRes = xxEnableRotorTempCoeffRes;

        end

        function [AmatND, BmatND, ResMatND, TnodesInit, ...
                  AdjacencyMat, InletCoolIdxs, OutletCoolIdxs] = generateSimulinkReducedOrderModel(obj, modelName, ...
                         coolingSystemsEnabled, BkptsStruct)
            % Generate a Simulink model implementing a ROM consisting of a 
            % set of state-space models at certain breakpoints (shaft speeds, 
            % coolant flow rates, and coolant inlet temperatures). The ROM
            % interpolates the state-space arrays between breakpoints.
            
            obj.turnOffLossDependenceWithTemperatureOrSpeed();

            % 1) COMPUTE ROM DATA
            
            % Enable the cooling systems and get the flowrate and inlet
            % temperature property names.
            disp("Enabling the specified cooling systems...");
            numCoolSys = length(coolingSystemsEnabled);
            frPropNames = cell(numCoolSys,1);
            TinPropNames = cell(numCoolSys,1);
            for idxCoolingSys = 1:length(string(coolingSystemsEnabled))
                thisCoolingSys = coolingSystemsEnabled{idxCoolingSys};
                enablePropName = strcat(erase(thisCoolingSys, ' '), '_Enable');
                obj.(enablePropName) = 1;

                frPropNames{idxCoolingSys} = strcat(erase(thisCoolingSys, ' '), '_FlowRate_m3ps');
                TinPropNames{idxCoolingSys} = strcat(erase(thisCoolingSys, ' '), '_InletTemperature_degC');            
            end
            obj.updateModel();

            % Get ambient temperature and initial node temperatures
            Tambient_degC =  obj.Tambient_degC;
            TnodesInit = Tambient_degC*ones(size(obj.TempMat));

            % size pre-allocation
            numStates = length(obj.TempMat);
            sizeMatND = nan(1, 1+2*numCoolSys+2);
            sizeMatND(1) = length(BkptsStruct.w);
            currIdx=2;
            for idxCool = 1:numCoolSys
                frNumel = length(BkptsStruct.(strcat('fr',num2str(idxCool))));
                sizeMatND(currIdx)=frNumel;
                currIdx=currIdx+1;
            end
            for idxCool = 1:numCoolSys
                TinNumel = length(BkptsStruct.(strcat('Tin',num2str(idxCool))));
                sizeMatND(currIdx)=TinNumel;
                currIdx=currIdx+1;
            end
            sizeMatND(end-1)=numStates;
            sizeMatND(end)=numStates;
            AmatND = nan(sizeMatND);
            BmatND = nan(sizeMatND);
            ResMatND = nan(sizeMatND);
            dummySS = ss(ones(numStates), ones(numStates), eye(numStates), zeros(numStates));
            StateSpaceND = repmat(dummySS, [1,1,sizeMatND(1:end-2)]);
            samplingGridStruct = getSamplingGridFromBkpts(BkptsStruct);
            StateSpaceND.SamplingGrid = samplingGridStruct;
            
            % Get coolant data
            AdjacencyMat = obj.AdjacencyMat;
            InletCoolIdxs = obj.InletArrayIdxs;
            OutletCoolIdxs = obj.OutletArrayIdxs;
            
            % Calculate state-space model at each breakpoint -------
            idxsCell = getNestedForLoopIdxs(sizeMatND(1:end-2));

            % Calculate state-space matrices for each breakpoint
            disp("Calculating state-space matrices for each breakpoint...");
            for idxBkptCmb = 1:length(idxsCell)
                bkptIdxsComb = idxsCell{idxBkptCmb};
                idxSpeed = bkptIdxsComb(1);
                disp("Breakpoint #" + num2str(idxBkptCmb));
                % Set speed
                disp("  Speed = " + num2str(BkptsStruct.w(idxSpeed)) + " rpm");
                obj.Shaft_Speed_RPM = BkptsStruct.w(idxSpeed);
                % Set cooling systems fr and Tin
                for idxCool = 1:numCoolSys
                    frBkpts = BkptsStruct.(strcat('fr',num2str(idxCool)));
                    TinBkpts = BkptsStruct.(strcat('Tin',num2str(idxCool)));
                    idxFr = idxCool+1;
                    idxTin = idxCool+numCoolSys+1;
                    frVal = frBkpts(bkptIdxsComb(idxFr));
                    disp("  " + coolingSystemsEnabled{idxCool} + " flow rate = " + frVal + " lpm");
                    obj.(frPropNames{idxCool}) = frVal/60/1000;
                    TinVal = TinBkpts(bkptIdxsComb(idxTin));
                    disp("  " + coolingSystemsEnabled{idxCool} +" inlet temperature = " + TinVal + " degC");
                    obj.(TinPropNames{idxCool}) = TinVal;
                end
                disp("");

                obj.updateModel();

                % get state-space matrices
                xAmat = obj.Amat;
                xBmat = obj.Bmat;
                xResMat = obj.ResMat;
                idxCell = num2cell(bkptIdxsComb);
                StateSpaceND(:,:,idxCell{:}) = ss(xAmat, xBmat, eye(numStates), zeros(numStates));
                for idx1 = 1:numStates
                    for idx2 = 1:numStates
                        indCell = num2cell([bkptIdxsComb,idx1,idx2],1);
                        ind = sub2ind(sizeMatND, indCell{:});
                        AmatND(ind) = xAmat(idx1,idx2);
                        BmatND(ind) = xBmat(idx1,idx2);
                        ResMatND(ind) = xResMat(idx1,idx2);
                    end
                end               
            end       

            % re-interpolate tables in a grid
            [ShaftTorqueVec, SpeedVec, xStator_Copper_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Stator_Copper_Loss_Mat);
            [~,~,xRotor_Cage_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Rotor_Cage_Loss_Mat);
            [~,~,xIron_Loss_Stator_Back_Iron_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Stator_Back_Iron_Mat);
            [~,~,xIron_Loss_Stator_Tooth_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Stator_Tooth_Mat);
            [~,~,xStray_Load_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Stray_Load_Loss_Mat);
            [~,~,xMagnet_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Magnet_Loss_Mat);
            [~,~,xIron_Loss_Rotor_Pole_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Rotor_Pole_Mat);
            [~,~,xIron_Loss_Rotor_Back_Iron_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Rotor_Back_Iron_Mat);
            [~,~,xIron_Loss_Rotor_Tooth_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Rotor_Tooth_Mat);
            [~,~,xFriction_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Friction_Loss_Mat);
            [~,~,xWindage_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Windage_Loss_Mat);  
            [~,~,xStator_Copper_Loss_AC_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Stator_Copper_Loss_AC_Mat);  
            [~,~,xBanding_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Banding_Loss_Mat);  
            [~,~,xSleeve_Loss_Mat] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Sleeve_Loss_Mat);  


            % Get PowerLossDistributor params
            disp("Estimating power loss distribution for each loss type...")
            lossDistrForEachType = obj.getLossDistrForEachType();
            TrefStator = obj.TrefStator_degC;
            TrefRotor = obj.TrefRotor_degC;
            xStatorCopperTempCoefResistivity = obj.StatorTempCoeffRes;
            xRotorCopperTempCoefResistivity = obj.RotorTempCoeffRes;
            xCoolantArrayIdxs = obj.CoolantArrayIdxs;            
            xCapMat = obj.CapMat;

            % 2) GENERATE SIMULINK MODEL
            disp("Creating Simulink model...")
            % create and open the model
            open_system(new_system(modelName));
            set_param(modelName, 'StopTime', '1000');
            % set solver params
            set_param(modelName,'AbsTol','1e-3'); 
            set_param(modelName,'AutoScaleAbsTol','off'); 
            % add blocks      
            [~,motName,~] = fileparts(obj.motFullFile);           
            romSubsysPath = strcat(modelName, '/',motName, '_ReducedOrderModel');
            add_block('built-in/Subsystem', romSubsysPath, ...
                'Position', '[0,0,250,150]');
            hLossMaps = add_block(strcat(obj.mcadROMLibName, '/LossMaps'), ...
                strcat(romSubsysPath, '/Loss Maps from Motor-CAD Lab'), ...
                'SpeedVec', 'SpeedVec', ...
                'ShaftTorqueVec', 'ShaftTorqueVec', ...
                'Stator_Copper_Loss_Mat', 'Stator_Copper_Loss_Mat', ...
                'Rotor_Cage_Loss_Mat', 'Rotor_Cage_Loss_Mat', ...
                'Iron_Loss_Stator_Back_Iron_Mat', 'Iron_Loss_Stator_Back_Iron_Mat', ...
                'Iron_Loss_Stator_Tooth_Mat', 'Iron_Loss_Stator_Tooth_Mat', ...
                'Stray_Loss_Stator_Iron_Proportion', 'Stray_Loss_Stator_Iron_Proportion', ...
                'Stray_Load_Loss_Mat', 'Stray_Load_Loss_Mat', ...
                'Magnet_Loss_Mat', 'Magnet_Loss_Mat', ...
                'Iron_Loss_Rotor_Pole_Mat', 'Iron_Loss_Rotor_Pole_Mat', ...
                'Iron_Loss_Rotor_Back_Iron_Mat', 'Iron_Loss_Rotor_Back_Iron_Mat', ...
                'Iron_Loss_Rotor_Tooth_Mat', 'Iron_Loss_Rotor_Tooth_Mat', ...
                'Friction_Loss_Mat', 'Friction_Loss_Mat', ...
                'Windage_Loss_Mat', 'Windage_Loss_Mat', ...
                'Stator_Copper_Loss_AC_Mat', 'Stator_Copper_Loss_AC_Mat', ...
                'Banding_Loss_Mat', 'Banding_Loss_Mat', ...
                'Sleeve_Loss_Mat', 'Sleeve_Loss_Mat' ...
                );
            hInTorque = add_block('simulink/Commonly Used Blocks/In1', strcat(romSubsysPath, '/ShaftTorque_Nm'));
            hInSpeed = add_block('simulink/Commonly Used Blocks/In1', strcat(romSubsysPath, '/ShaftSpeed_RPM')); 
            hOutTnodes = add_block('simulink/Commonly Used Blocks/Out1', strcat(romSubsysPath, '/TempNodes'));
            if isempty(coolingSystemsEnabled) % passive cooling
                hIntpSs = add_block(strcat(obj.mcadROMLibName, '/InterpolatedStateSpaceThermalModel (LPV) (Passive Cooling)'), ...
                    strcat(romSubsysPath, '/nterpolatedStateSpaceThermalModel (LPV) (Passive Cooling)'), ...
                    'StateSpaceND', 'StateSpaceND', ...
                    'TnodesInit', 'TnodesInit'...
                    );
                hPowLossDistr = add_block(strcat(obj.mcadROMLibName, '/PowerLossDistributor (Passive Cooling)'), ...
                    strcat(romSubsysPath, '/PowerLossDistributor (Passive Cooling)'), ...
                    'TrefStator', 'TrefStator', ...
                    'StatorCopperTempCoefResistivity', 'StatorCopperTempCoefResistivity', ...
                    'TrefRotor', 'TrefRotor', ...
                    'RotorCopperTempCoefResistivity', 'RotorCopperTempCoefResistivity', ...
                    'LossDistrForEachType', 'LossDistrForEachType', ...
                    'CapMat', 'CapMat');                
            else % active cooling
                hIntpSs = add_block(strcat(obj.mcadROMLibName, '/InterpolatedStateSpaceThermalModel (LPV)'), ...
                    strcat(romSubsysPath, '/Interpolated State-Space Thermal Model (LPV)'), ...
                    'StateSpaceND', 'StateSpaceND', ...
                    'CoolAdjMat', 'AdjacencyMat', ...
                    'CoolantArrayIdxs', 'CoolantArrayIdxs', ...
                    'InletCoolantIdxs', 'InletCoolIdxs', ...
                    'OutletCoolantIdxs', 'OutletCoolIdxs', ...
                    'TnodesInit', 'TnodesInit'...
                    );
                hPowLossDistr = add_block(strcat(obj.mcadROMLibName, '/PowerLossDistributor'), ...
                    strcat(romSubsysPath, '/Power Loss Distributor'), ...
                    'TrefStator', 'TrefStator', ...
                    'StatorCopperTempCoefResistivity', 'StatorCopperTempCoefResistivity', ...
                    'TrefRotor', 'TrefRotor', ...
                    'RotorCopperTempCoefResistivity', 'RotorCopperTempCoefResistivity', ...
                    'LossDistrForEachType', 'LossDistrForEachType', ...
                    'CapMat', 'CapMat', ...
                    'CoolantArrayIdxs', 'CoolantArrayIdxs');
                hDemuxMux = add_block(strcat(obj.mcadROMLibName, '/CoolantDemuxMux'), ...
                    strcat(romSubsysPath, '/Coolant Interface'), ...
                    'numCoolSys', num2str(numCoolSys));
                hCoolIn = [];          
                for idxCool = 1:numCoolSys
                    thisCoolingSys = coolingSystemsEnabled{idxCool};
                    inletName = strcat(erase(thisCoolingSys, ' '), '_Inlet');
                    hCoolIn(end+1) = add_block('simulink/Commonly Used Blocks/In1', strcat(romSubsysPath, '/', inletName)); %#ok<AGROW> 
                end
                hOutPnodes = add_block('simulink/Commonly Used Blocks/Out1', strcat(romSubsysPath, '/PowNodes'));
                hOutTcoolOut = add_block('simulink/Commonly Used Blocks/Out1', strcat(romSubsysPath, '/TcoolOutVec'));
            end

            % connect blocks
            hLossMapsPorts = get_param(hLossMaps, 'PortHandles');
            hInTorquePort = get_param(hInTorque, 'PortHandles');
            hInSpeedPort = get_param(hInSpeed, 'PortHandles');
            hPowLossDistrPorts = get_param(hPowLossDistr, 'PortHandles');
            hIntpSsPorts = get_param(hIntpSs, 'PortHandles');
            hOutTnodesPort = get_param(hOutTnodes, 'PortHandles');
            add_line(romSubsysPath, hInTorquePort.Outport, hLossMapsPorts.Inport(1));
            add_line(romSubsysPath, hInSpeedPort.Outport, hLossMapsPorts.Inport(2)); 
            add_line(romSubsysPath, hPowLossDistrPorts.Outport, hIntpSsPorts.Inport(1));
            add_line(romSubsysPath, hInSpeedPort.Outport, hIntpSsPorts.Inport(2));
            add_line(romSubsysPath, hIntpSsPorts.Outport(1), hPowLossDistrPorts.Inport(1));
            add_line(romSubsysPath, hLossMapsPorts.Outport, hPowLossDistrPorts.Inport(2));
            add_line(romSubsysPath, hIntpSsPorts.Outport(1), hOutTnodesPort.Inport);
            if ~isempty(coolingSystemsEnabled) % active cooling  - specific additional connections
                hDemuxMuxPorts = get_param(hDemuxMux, 'PortHandles');
                hOutPnodesPort = get_param(hOutPnodes, 'PortHandles');
                hOutTcoolOutPort = get_param(hOutTcoolOut, 'PortHandles');         
                for idxCool = 1:numCoolSys
                    hCoolInPort = get_param(hCoolIn(idxCool), 'PortHandles');    
                    add_line(romSubsysPath, hCoolInPort.Outport, hDemuxMuxPorts.Inport(idxCool));
                end
                add_line(romSubsysPath, hDemuxMuxPorts.Outport(1), hIntpSsPorts.Inport(3));
                add_line(romSubsysPath, hDemuxMuxPorts.Outport(2), hIntpSsPorts.Inport(4));
                add_line(romSubsysPath, hIntpSsPorts.Outport(2), hOutPnodesPort.Inport);
                add_line(romSubsysPath, hIntpSsPorts.Outport(3), hOutTcoolOutPort.Inport);
            end
            
            % Write data into model workspace
            mdlWks = get_param(modelName,'ModelWorkspace');
            assignin(mdlWks,'StateSpaceND', StateSpaceND);
            assignin(mdlWks,'AdjacencyMat', AdjacencyMat);
            assignin(mdlWks,'InletCoolIdxs', InletCoolIdxs);
            assignin(mdlWks,'OutletCoolIdxs', OutletCoolIdxs);
            assignin(mdlWks,'TnodesInit', TnodesInit);
            assignin(mdlWks,'SpeedVec', SpeedVec);
            assignin(mdlWks,'ShaftTorqueVec', ShaftTorqueVec);
            assignin(mdlWks,'Stator_Copper_Loss_Mat', xStator_Copper_Loss_Mat);
            assignin(mdlWks,'Rotor_Cage_Loss_Mat', xRotor_Cage_Loss_Mat);
            assignin(mdlWks,'Iron_Loss_Stator_Back_Iron_Mat', xIron_Loss_Stator_Back_Iron_Mat);
            assignin(mdlWks,'Iron_Loss_Stator_Tooth_Mat', xIron_Loss_Stator_Tooth_Mat); 
            assignin(mdlWks,'Stray_Loss_Stator_Iron_Proportion', obj.Stray_Loss_Stator_Iron_Proportion); 
            assignin(mdlWks,'Stray_Load_Loss_Mat', xStray_Load_Loss_Mat); 
            assignin(mdlWks,'Magnet_Loss_Mat', xMagnet_Loss_Mat);
            assignin(mdlWks,'Iron_Loss_Rotor_Pole_Mat', xIron_Loss_Rotor_Pole_Mat);
            assignin(mdlWks,'Iron_Loss_Rotor_Back_Iron_Mat', xIron_Loss_Rotor_Back_Iron_Mat);
            assignin(mdlWks,'Iron_Loss_Rotor_Tooth_Mat', xIron_Loss_Rotor_Tooth_Mat);
            assignin(mdlWks,'Friction_Loss_Mat', xFriction_Loss_Mat);
            assignin(mdlWks,'Windage_Loss_Mat', xWindage_Loss_Mat);
            assignin(mdlWks,'Stator_Copper_Loss_AC_Mat', xStator_Copper_Loss_AC_Mat);
            assignin(mdlWks,'Banding_Loss_Mat', xBanding_Loss_Mat);
            assignin(mdlWks,'Sleeve_Loss_Mat', xSleeve_Loss_Mat);
            assignin(mdlWks,'TrefStator', TrefStator);
            assignin(mdlWks,'TrefRotor', TrefRotor);
            assignin(mdlWks,'StatorCopperTempCoefResistivity', xStatorCopperTempCoefResistivity);
            assignin(mdlWks,'RotorCopperTempCoefResistivity', xRotorCopperTempCoefResistivity);
            assignin(mdlWks,'LossDistrForEachType', lossDistrForEachType);
            assignin(mdlWks,'CapMat', xCapMat);
            assignin(mdlWks,'CoolantArrayIdxs', xCoolantArrayIdxs);

            % add dummy inputs to ROM subsystem
            add_block('simulink/Sources/Constant', strcat(modelName, '/SpeedRPM'), ...
                'Value', '2000');
            add_block('simulink/Sources/Constant', strcat(modelName, '/TorqueNm'), ...
                'Value', '20');
            add_line(modelName, strcat('TorqueNm', '/1'), strcat(motName, '_ReducedOrderModel', '/1'))
            add_line(modelName, strcat('SpeedRPM', '/1'), strcat(motName, '_ReducedOrderModel', '/2'))
            for idxCool = 1:numCoolSys
                thisCoolingSys = coolingSystemsEnabled{idxCool};
                frBlockName = strcat(erase(thisCoolingSys, ' '), '_Flowrate_lpm');
                TinBlockName = strcat(erase(thisCoolingSys, ' '), '_InletTemp_degC');
                add_block('simulink/Sources/Constant', strcat(modelName, '/', frBlockName), ...
                'Value', '3');
                add_block('simulink/Sources/Constant', strcat(modelName, '/', TinBlockName), ...
                'Value', '25');
                add_block('simulink/Commonly Used Blocks/Mux', strcat(modelName, '/Mux', num2str(idxCool)));
                add_line(modelName, strcat(frBlockName, '/1'), strcat('Mux', num2str(idxCool), '/1'))
                add_line(modelName, strcat(TinBlockName, '/1'), strcat('Mux', num2str(idxCool), '/2'))
                add_line(modelName, strcat('Mux', num2str(idxCool), '/1'), strcat(motName, '_ReducedOrderModel', '/', num2str(2+idxCool)));               
            end

            % add scopes and outport to the ROM subsystem outputs
            add_block('simulink/Sinks/Scope', strcat(modelName, '/Node Temperatures'))
            add_block('simulink/Sinks/Out1', strcat(modelName, '/Out'))
            add_line(modelName, strcat(motName, '_ReducedOrderModel', '/1'), strcat('Node Temperatures', '/1'))
            add_line(modelName, strcat(motName, '_ReducedOrderModel', '/1'), strcat('Out', '/1'))
            if ~isempty(coolingSystemsEnabled) % active cooling  - specific additional Scopes
                add_block('simulink/Sinks/Scope', strcat(modelName, '/Node Powers'))
                add_block('simulink/Sinks/Scope', strcat(modelName, '/Coolant Outlet Temperature'))
                add_line(modelName, strcat(motName, '_ReducedOrderModel', '/2'), strcat('Node Powers', '/1'))
                add_line(modelName, strcat(motName, '_ReducedOrderModel', '/3'), strcat('Coolant Outlet Temperature', '/1'))
            end
            
            % tidy-up connections
            Simulink.BlockDiagram.arrangeSystem(romSubsysPath);
            Simulink.BlockDiagram.arrangeSystem(modelName);          

        end

    end

    methods(Access=private)
        function setupForThermalMatricesUpdate(obj)
            % Matrices must be computed at zero loss inputs and no
            % temperature coefficient effect
            
            obj.EnableStatorTempCoeffRes = 0; % no temp coefficient effect
            obj.EnableRotorTempCoeffRes = 0; % no temp coefficient effect
            residualLossVal = 0.1; % very small loss ~ 0
            lossVec = residualLossVal*ones(size(obj.LossValues)); 
            
            obj.runThermalSteadyStateWithSpecifiedLosses(lossVec)

        end

        function turnOffLossDependenceWithTemperatureOrSpeed(obj)
            % Disable loss dependence on speed and temperature. This is
            % required to correctly estimate the loss distribution amongst
            % nodes.
            if obj.Loss_Function_Speed==int32(1)
                warning("Turning off 'Copper Losses Vary with Temperature'. This is required to correctly estimate the loss distribution amongst nodes.")
                obj.Loss_Function_Speed = 0;
            end

            if obj.Copper_Losses_Vary_With_Temperature==int32(1)
                warning("Turning off 'Copper Losses Vary with Temperature'. This is required to correctly estimate the loss distribution amongst nodes.")
                obj.Copper_Losses_Vary_With_Temperature = 0;
            end

            if obj.RotorCopperLossesVaryWithTemp==int32(1)
                warning("Turning off 'Rotor Cage Losses Vary with Temperature'. This is required to correctly estimate the loss distribution amongst nodes.")
                obj.RotorCopperLossesVaryWithTemp = 0;
            end

            if obj.StatorIronStrayLoadLossesVaryWithTemp==int32(1)
                warning("Turning off 'Stator Iron Stray Losses Vary with Temperature'. This is required to correctly estimate the loss distribution amongst nodes.")
                obj.StatorIronStrayLoadLossesVaryWithTemp = 0;
            end

            if obj.RotorIronStrayLoadLossesVaryWithTemp==int32(1)
                warning("Turning off 'Rotor Iron Stray Losses Vary with Temperature'. This is required to correctly estimate the loss distribution amongst nodes.")
                obj.RotorIronStrayLoadLossesVaryWithTemp = 0;
            end

            if obj.StatorCopperStrayLoadLossesVaryWithTemp==int32(1)
                warning("Turning off 'Stator Copper Stray Losses Vary with Temperature'. This is required to correctly estimate the loss distribution amongst nodes.")
                obj.StatorCopperStrayLoadLossesVaryWithTemp = 0;
            end

            if obj.RotorCopperStrayLoadLossesVaryWithTemp==int32(1)
                warning("Turning off 'Rotor Copper Stray Losses Vary with Temperature'. This is required to correctly estimate the loss distribution amongst nodes.")
                obj.RotorCopperStrayLoadLossesVaryWithTemp = 0;
            end

        end
    end

end

% Helper functions

function [Amat, Bmat] = getStateSpaceMatricesFromThermalMatrices(CapMat, ResMat)
    % Get A,B matrices from capacitance and resistance matrices

    %  Create matrixes A B C and D for state space model
    [N,~]=size(ResMat);
    R_inv = zeros(N);
    C_inv = zeros(N,1);
    if numel(CapMat)>N % CapMat is a NxN matrix
        CapVec = diag(CapMat);
        CapVec = CapVec(:); % column vec
    else
        CapVec = CapMat(:);
    end
    
    % Setup model parameters
    % Build inverse R,C matrixes
    for row=1:1:N
        for col=1:1:N
            R_inv(row,col)=1/ResMat(row,col);
            if row==col
                R_inv(row,col)=0;
            end
        end
    end
    
    for i=1:1:N
        C_inv(i,1)=1/CapVec(i,1);
    end
    
    % Build Amat matrix
    
    Amat=zeros(N);
    
    % Top right half of Amat matrix
    for row=1:1:N
        for col=row+1:1:N
          if C_inv(row,1) < 10^6           
            Amat(row,col)=C_inv(row,1)*R_inv(row,col);
          else
            Amat(row,col)=(10^6)*R_inv(row,col);
          end
        end
    end
    
    % Bottom left half of Amat matrix
    for row=2:1:N
        for col=1:row-1    
            if C_inv(row,1)< 10^6
                Amat(row,col)=C_inv(row,1)*R_inv(row,col);
            else
                Amat(row,col)=(10^6)*R_inv(row,col);
            end
        end
    end
    
    % Diagonal of Amat Matrix
    for row=1:1:N 
      if C_inv(row,1)< 10^6
        Amat(row,row)=-(C_inv(row,1)*(sum(R_inv(row,:))));
      else
        Amat(row,row)=(-(10^6)*(sum(R_inv(row,:))));
      end
    end
    
    % Build Bmat matrix
    Bmat=zeros(N);
    
    for row=1:1:N  
        if C_inv(row,1)< 10^6
            Bmat(row,row)=C_inv(row,1);   
    
        else 
            Bmat(row,row)= 10^6;
        end   
    end

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
    % Read the resistance values. Extra node for ambient.
    for i=1:numNodes
        formatSpec = [formatSpec, '%f']; %#ok<AGROW> 
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

function idxsCell = getNestedForLoopIdxs(arraySizeVec)

    n = length(arraySizeVec);   % Number of indices
    currentIdxs = ones(1, n);   % Current indices
    idxsCell = cell(prod(arraySizeVec), 1);
    counter = 0;
    ready = false;
    while ~ready
        counter = counter + 1;
        idxsCell{counter} = currentIdxs;
        % Increase the index vector
        ready = true;    
        for k = n:-1:1
          currentIdxs(k) = currentIdxs(k) + 1;
          if currentIdxs(k) <= arraySizeVec(k)
            ready = false;
            break;
          end
          currentIdxs(k) = 1;  % Reset, proceed with previous element
        end
    end

end

function [xVec, yVec, zgMat] = reInterpolateTable(xsMat, ysMat, zsMat)
        % Re-interpolate scattered map into square grid

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

function samplingGridStruct = getSamplingGridFromBkpts(BkptsStruct)
    % Get sampling grid from breakpoints structure.

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
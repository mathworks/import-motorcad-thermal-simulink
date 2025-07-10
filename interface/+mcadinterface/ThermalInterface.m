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

        CapMat (:,1) double % Node capacitance vector
        ResMat (:,:) double % Node resistance matrix
        PowMat (:,1) double % Node steady-state power vector
        TempMat (:,1) double % Node temperature boundary condition

        % State-space matrices

        Amat (:,:) double % State-space A matrix
        Bmat (:,:) double % State-space B matrix

        % Cooling system data      

        EnabledCoolingSystems (:,1) logical % List of which cooling systems are enabled (true or false)
        CoolingSystemsDigraphs (:,1) struct % Structure containing the directed graph of each cooling systems node flow
        AdjacencyMat (:,:) double % Global adjacency matrix for cooling systems node flow connectivity
        InletArrayIdxs (:,1) cell % Inlet nodes indices
        OutletArrayIdxs (:,1) cell % Outlet nodes indices
        CoolantArrayIdxs (:,1) double % Coolant nodes indices
    end

    properties(Constant)
        % List of supported cooling systems
        SupportedCoolingSystemsMcadNames = ...
                    {'Blown Over'; ...
                     'Ventilated'; ...
                     'Housing Water Jacket'; ...
                     'Shaft Spiral Groove'; ...
                     'Wet Rotor'; ...
                     'Spray Cooling'; ... % standard spray cooling
                     'Spray_RadialHousing_F'; ... % multi-nozzle spray - radial hosing (front side) 
                     'Spray_RadialHousing_R'; ... % multi-nozzle spray - radial hosing (rear side) 
                     'Spray_RadialRotor_F'; ...   % multi-nozzle spray - radial rotor (front side) 
                     'Spray_RadialRotor_R'; ...   % multi-nozzle spray - radial rotor (rear side) 
                     'Spray_AxialEndcap_F'; ...   % multi-nozzle spray - axial endcap (front side) 
                     'Spray_AxialEndcap_R'; ...   % multi-nozzle spray - axial endcap (rear side) 
                     'Rotor Water Jacket'; ...
                     'Slot Water Jacket'; ...
                                   }; 
    end

    properties(Constant, Access=private)
        MultiNozzleSprayCoolingSystems=["Spray_RadialHousing_F"; ...
                                        "Spray_RadialHousing_R"; ...
                                        "Spray_RadialRotor_F"; ...
                                        "Spray_RadialRotor_R"; ...
                                        "Spray_AxialEndcap_F"; ...
                                        "Spray_AxialEndcap_R"; ...
                                        ];
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
            % obj.checkCoolingSystemsIndependent();

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

        function [CapMat, ResMat, PowMat, TempMat, McadIdxes] = getThermalStateSpaceMatricesFromFiles(obj)
            %GETTHERMALSTATEMATRIXSFROMFILES Retrieves thermal matrices and corresponding Motor-CAD indices from Motor-CAD export files.
            %
            % [CapMat, ResMat, PowMat, TempMat, McadIdxes] = getThermalStateSpaceMatricesFromFiles(obj)
            %
            % Outputs:
            %   CapMat    - (Nx1 double) Node capacitance vector.
            %   ResMat    - (NxN double) Node resistance matrix.
            %   PowMat    - (Nx1 double) Node steady-state power vector.
            %   TempMat   - (Nx1 double) Node temperature boundary condition.
            %   McadIdxes - (Nx1 int32) List of corresponding Motor-CAD indices.
        
            motFileName = obj.motFullFile;
        
            CapFileName = strrep(motFileName, '.mot', '.cmf');
            ResFileName = strrep(motFileName, '.mot', '.rmf');
            PowFileName = strrep(motFileName, '.mot', '.pmf');
            TempFileName = strrep(motFileName, '.mot', '.tmf');
        
            % Read cmf file
            [CapMatNm1, CapMcadIdxes] = readMfFile1D(CapFileName); % size numNodes-1
            numNodes = length(CapMatNm1) + 1; % number of nodes (+1 for ambient)
            CapMat = zeros(numNodes,1);
            CapMat(2:end) = CapMatNm1;
        
            % Read rmf file
            ResMat = readMfFile2D(ResFileName, numNodes); % Assuming readMfFile2D remains unchanged
        
            % Read pmf file
            [PowMatNm1, PowMcadIdxes] = readMfFile1D(PowFileName); % size numNodes-1
            PowMat = zeros(numNodes,1);
            PowMat(2:end) = PowMatNm1;
        
            % Read tmf file
            [TempMat, TempMcadIdxes] = readMfFile1D(TempFileName); % size numNodes
        
            % Handle Ambient Node
            CapMat(1,1) = 1e20; % Ambient must have large capacitance
            for idx_node = 2:numNodes      
                if TempMat(idx_node) ~= -10000000 % Fixed-temperature node
                    % Model these nodes as having very large capacitances
                    CapMat(idx_node) = 1e20;
                end 
            end
        
            % Consistency Check: Ensure that Cap, Pow, and Temp files have consistent indices
            % Assuming that Motor-CAD indices are unique and ordered, except for ambient
            % If not, additional mapping may be required
        
            % For simplicity, assuming that CapMcadIdxes and PowMcadIdxes correspond to nodes 2:numNodes
            % and TempMcadIdxes includes node 1 (ambient) followed by nodes 2:numNodes
        
            % Assign Motor-CAD indices
            % Ambient node assumed to have index 0 or a special identifier
            % Here, we will assign a special index for ambient, e.g., 0
            McadIdxes = zeros(numNodes,1); % Initialize with 0 for ambient
            McadIdxes(2:end) = CapMcadIdxes; % Assign indices from .cmf (excluding ambient)
        
            % Optionally, verify that PowMcadIdxes match CapMcadIdxes
            if ~isequal(PowMcadIdxes, CapMcadIdxes)
                warning('Mismatch between Motor-CAD indices in .cmf and .pmf files.');
            end
            if ~isequal(TempMcadIdxes(2:end), CapMcadIdxes)
                warning('Mismatch between Motor-CAD indices in .cmf and .tmf files.');
            end
        
        end

        function [GroupNamesAndMcadIdxes, NodeNamesAndMcadIdx] = updateMatricesAndGroupNamesAndNodes(obj)
            %UPDATEMATRICESANDGROUPNAMESANDNODES Updates state-space matrices, group names, and node data 
            %based on current Motor-CAD state.
            %
            % [GroupNamesAndMcadIdxes, NodeNamesAndMcadIdx] = updateMatricesAndGroupNamesAndNodes(obj)
            %
            % Outputs:
            %   GroupNamesAndMcadIdxes - Cell array with group names and their associated Motor-CAD node indices.
            %   NodeNamesAndMcadIdx    - Cell array with node names and their corresponding Motor-CAD node indices.
        
            % --- 1) BACK UP CURRENT PROPERTIES (SO YOU CAN RESTORE AFTER MATRIX EXPORT) ---
            xxLossValues = obj.LossValues;
            xxEnableStatorTempCoeffRes = obj.EnableStatorTempCoeffRes;
            xxEnableRotorTempCoeffRes = obj.EnableRotorTempCoeffRes;
        
            % --- 2) GENERATE (TEMPORARY) THERMAL MATRICES & WRITE FILES ---
            obj.setupForThermalMatricesUpdate();
            obj.writeThermalStateSpaceFiles();
            [xCapMat, xResMat, xPowMat, xTempMat, xMcadIdxs] = obj.getThermalStateSpaceMatricesFromFiles();
        
            % --- 3) RESTORE PROPERTIES AFTER MATRIX EXPORT ---
            obj.LossValues = xxLossValues; 
            obj.EnableStatorTempCoeffRes = xxEnableStatorTempCoeffRes;
            obj.EnableRotorTempCoeffRes = xxEnableRotorTempCoeffRes;
        
            % --- 4) UPDATE THERMAL MATRICES & STATE-SPACE ARRAYS IN THE OBJECT ---
            obj.CapMat = xCapMat;
            obj.ResMat = xResMat;
            obj.PowMat = xPowMat;
            obj.TempMat = xTempMat;
            [xAmat, xBmat] = getStateSpaceMatricesFromThermalMatrices(xCapMat, xResMat);
            obj.Amat = xAmat;
            obj.Bmat = xBmat;
        
            % --- 5) READ THE .NMF FILE (WHERE NODE NAMES ARE NOT TRUNCATED) ---
            motFileName = obj.motFullFile;
            NodeFileName = strrep(motFileName, '.mot', '.nmf');
            [GroupNames, nodeNames, GroupIdxs, nmfMcadIdxs] = readNmfFile(NodeFileName);
            %
            %  > GroupNames  : e.g. {'Stator (Active)', 'Spray Cooling', ...}
            %  > nodeNames   : e.g. {'StatorSlot1', 'Spray_RadialHousing_Inlet_Fluid_F', ...}
            %  > GroupIdxs   : each cell is a numeric array of MCAD node indices in that group
            %  > nmfMcadIdxs : cell array of MCAD node indices (one for each entry in nodeNames)
        
            % --- 6) BUILD GROUP->IDX MAPPING ---
            numGroups = length(GroupNames);
            GroupNamesAndMcadIdxes = cell(numGroups, 2);
            for g = 1:numGroups
                GroupNamesAndMcadIdxes{g, 1} = GroupNames{g};
                GroupNamesAndMcadIdxes{g, 2} = GroupIdxs{g}; 
            end
        
            % --- 7) BUILD NODE->IDX MAPPING, THEN RE-SORT TO MATCH xMcadIdxs ORDER ---
            numNodes = length(nmfMcadIdxs);
            NodeNamesAndMcadIdx = cell(numNodes, 2);
        
            % Populate initial (unsorted) mapping from .nmf
            for n = 1:numNodes
                NodeNamesAndMcadIdx{n,1} = nodeNames{n};      % Full name
                NodeNamesAndMcadIdx{n,2} = nmfMcadIdxs{n};    % Numeric index
            end
        
            % We want to reorder so that NodeNamesAndMcadIdx(:) align with the row/col
            % order in xCapMat, xResMat, xPowMat, xTempMat, i.e. the order of xMcadIdxs.
            sortedNodeIdx = nan(numNodes,1);
            for i = 1:numNodes
                thisMcadIdx = xMcadIdxs(i);          % MCAD index from the matrix file
                % find where that index appears in nmfMcadIdxs
                pos = find(cell2mat(nmfMcadIdxs) == thisMcadIdx, 1);
                if isempty(pos)
                    error('Motor-CAD index %d from .cmf/.tmf/.pmf/.rmf not found in .nmf.', thisMcadIdx);
                end
                sortedNodeIdx(i) = pos;
            end
        
            % Reorder NodeNamesAndMcadIdx to match xMcadIdxs
            NodeNamesAndMcadIdx = NodeNamesAndMcadIdx(sortedNodeIdx, :);
        
            % --- 8) SET obj.NodeNames USING THE FULL (UNCROPPED) NAMES FROM .NMF ---
            obj.NodeNames = cell(numNodes,1);
            for i = 1:numNodes
                obj.NodeNames{i} = stripOuterParentheses(NodeNamesAndMcadIdx{i,1});
            end
        
            % --- 9) SAVE FINAL MAPPINGS BACK TO THE OBJECT ---
            obj.GroupNamesAndMcadIdxes = GroupNamesAndMcadIdxes;
            obj.NodeNamesAndMcadIdx    = NodeNamesAndMcadIdx;

        end

        function updateCoolingSystemsData(obj)
            % Refresh the adjacency, inlet, and outlet data for each cooling system.
        
            % 1) Split out the spray-cooling groups as needed:
            postProcessedGroups = splitSprayCoolingGroups(obj);

            % 2) Add the Blown-Over cooling
            postProcessedGroups{end+1,1} = 'Blown Over';
            postProcessedGroups{end,2} = [0, 0]; % inlet=ambient, outlet=ambient. Mcad ambient index is always 0
            
            % 3) Now filter them by the supported names:
            groupNames = postProcessedGroups(:,1);
            isSupported = contains(groupNames, obj.SupportedCoolingSystemsMcadNames);
            
            obj.CoolingSystemNamesAndMcadIdxes = postProcessedGroups(isSupported,:);
            
            % 4) Figure out which systems are actually "enabled" in the final object.
            obj.EnabledCoolingSystems = contains(obj.SupportedCoolingSystemsMcadNames, ...
                                                 obj.CoolingSystemNamesAndMcadIdxes(:,1));
            
            numSuppCoolSys = length(obj.SupportedCoolingSystemsMcadNames);
            obj.CoolingSystemsDigraphs = struct();
            
            % 5) Build each coolant's directed graph
            for idxSuppCoolSys = 1:numSuppCoolSys
                isEnabled = obj.EnabledCoolingSystems(idxSuppCoolSys);
                coolantName = obj.SupportedCoolingSystemsMcadNames{idxSuppCoolSys};
        
                if isEnabled && not(strcmp(coolantName, 'Blown Over'))
                    [digG, nodeInfo] = obj.getDigraphForCoolantGroup(coolantName, false);
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).Digraph              = digG;
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).NodeNamesAndArrayIdx = nodeInfo;
                else
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).Digraph              = [];
                    obj.CoolingSystemsDigraphs(idxSuppCoolSys).NodeNamesAndArrayIdx = [];
                end
            end
        
            % 6) Finally assemble global adjacency plus inlet/outlet indices
            [A, inletsCell, outletsCell, coolantIdxs] = getAdjacencyMatAndInletOutletIdxs(obj);
            obj.AdjacencyMat     = A;
            obj.InletArrayIdxs   = inletsCell;    % cell array, one cell per enabled system
            obj.OutletArrayIdxs  = outletsCell;   % likewise
            obj.CoolantArrayIdxs = coolantIdxs;
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
            %----------------------------------------------------------------------
            % 1) Identify all nodes belonging to this coolant group
            %----------------------------------------------------------------------
            AllNodeNames = obj.NodeNames;
            AllMcadIdxes = [obj.NodeNamesAndMcadIdx{:,2}];
            xResMat = obj.ResMat;
        
            % Find the actual MCAD indices for this group
            matchRow = strcmp(obj.CoolingSystemNamesAndMcadIdxes(:,1), CoolantGroupName);
            if ~any(matchRow)
                error('Cooling group "%s" not found among obj.CoolingSystemNamesAndMcadIdxes.', CoolantGroupName);
            end
            CoolantNodesMcadIdx = obj.CoolingSystemNamesAndMcadIdxes{matchRow,2};
        
            % Map MCAD indices to array indices
            numCoolantNodes = length(CoolantNodesMcadIdx);
            CoolantNodeNames = cell(numCoolantNodes,1);
            CoolantArrayIdxes = zeros(numCoolantNodes,1);
            NodeNamesAndArrayIdx = cell(numCoolantNodes,3);  % {NodeName, ArrayIdx, McadIdx}
        
            for iC = 1:numCoolantNodes
                thisMcadIdx = CoolantNodesMcadIdx(iC);
                arrIdx = find(AllMcadIdxes == thisMcadIdx);
                CoolantArrayIdxes(iC) = arrIdx;
                CoolantNodeNames{iC} = AllNodeNames{arrIdx};
        
                NodeNamesAndArrayIdx{iC,1} = CoolantNodeNames{iC};
                NodeNamesAndArrayIdx{iC,2} = arrIdx;
                NodeNamesAndArrayIdx{iC,3} = thisMcadIdx;
            end
        
            %----------------------------------------------------------------------
            % 2) Build the sub-matrix of resistances and create an undirected adjacency
            %----------------------------------------------------------------------
            ResMatCoolant = xResMat(CoolantArrayIdxes, CoolantArrayIdxes);
        
            % A small threshold to detect "connected" elements:
            UndirectedAdjacencyMatrix = (abs(ResMatCoolant) < 1e8) & (ResMatCoolant ~= 0);
        
            % Make sure it's symmetric, since we consider it an "undirected" view
            UndirectedAdjacencyMatrix = UndirectedAdjacencyMatrix | UndirectedAdjacencyMatrix';
        
            % Create the undirected graph from the sub-adjacency
            % (MATLAB requires it be symmetric, which we ensured above)
            UndirectedGraphCoolant = graph(UndirectedAdjacencyMatrix, CoolantNodeNames);
        
            %----------------------------------------------------------------------
            % 3) Identify all inlets by name
            %----------------------------------------------------------------------
            isInlet = contains(CoolantNodeNames, 'Inlet', 'IgnoreCase', true);
            inletNames = CoolantNodeNames(isInlet);
            if isempty(inletNames)
                warning('No node name contains "Inlet" in coolant group "%s". Assuming no inlets.', CoolantGroupName);
            end
        
            %----------------------------------------------------------------------
            % 4) Construct a directed graph with edges oriented away from each inlet
            %----------------------------------------------------------------------
            % Start with a fully bidirectional digraph:
            DirectedGraphCoolant = digraph(adjacency(UndirectedGraphCoolant), CoolantNodeNames);
        
            % For each inlet, do a BFS in the *undirected* graph to discover a tree,
            % then remove the back edge (child->parent) from the directed graph.
            for iInlet = 1:length(inletNames)
                thisInletName = inletNames{iInlet};
        
                % BFS on the undirected graph
                bfsNodeList = UndirectedGraphCoolant.bfsearch(thisInletName);
        
                % Consecutive pairs (parent->child) in BFSNodeList define the BFS tree edges.
                % We keep the direction parent->child, but remove child->parent.
                for k = 2:length(bfsNodeList)
                    parentNode = bfsNodeList(k-1);
                    childNode  = bfsNodeList(k);
        
                    if DirectedGraphCoolant.findedge(childNode, parentNode) > 0
                        DirectedGraphCoolant = DirectedGraphCoolant.rmedge(childNode, parentNode);
                    end
                end
            end
        
            %----------------------------------------------------------------------
            % 5) Identify potential outlets: any node with outdegree=0
            %    But we may discover multiple or none. We'll handle the "none" case below.
            %----------------------------------------------------------------------
            outDeg = outdegree(DirectedGraphCoolant);
            isOutlet = (outDeg == 0);
        
            candidateOutletNames = DirectedGraphCoolant.Nodes.Name(isOutlet);
        
            %----------------------------------------------------------------------
            % 6) If we discover *no* outlets, add a virtual outlet
            %    (applies, e.g., for pure spray cooling or other types with an unmodeled sink)
            %----------------------------------------------------------------------
            if isempty(candidateOutletNames)
                % Create a new node in the directed graph
                vName = 'VirtualOutlet';
                DirectedGraphCoolant = addnode(DirectedGraphCoolant, vName);
        
                % For every inlet, connect inlet -> VirtualOutlet
                for iInlet = 1:length(inletNames)
                    inName = inletNames{iInlet};
                    DirectedGraphCoolant = addedge(DirectedGraphCoolant, inName, vName, 1);
                end
        
                % Also update the NodeNamesAndArrayIdx to reflect the new node
                % We'll append an entry with a new "array index" and "mcad index"
                newArrayIdx = max(cell2mat(NodeNamesAndArrayIdx(:,2))) + 1;
                newMcadIdx  = max(cell2mat(NodeNamesAndArrayIdx(:,3))) + 1;
        
                NodeNamesAndArrayIdx{end+1,1} = vName;      % node name
                NodeNamesAndArrayIdx{end,2}   = newArrayIdx; 
                NodeNamesAndArrayIdx{end,3}   = newMcadIdx; 
        
                % No need to alter ResMatCoolant inside this function unless you 
                % want to track it in your overall adjacency (that can happen later).
            end
        
            %----------------------------------------------------------------------
            % 7) Remove cycles discovered by a DFS across all inlets
            %----------------------------------------------------------------------
            if ~isempty(inletNames)
                [~, eidx] = DirectedGraphCoolant.dfsearch(inletNames{1}, 'edgetodiscovered', 'Restart', true);
                DirectedGraphCoolant = DirectedGraphCoolant.rmedge(eidx);
            end
        
            %----------------------------------------------------------------------
            % 8) Plot if needed
            %----------------------------------------------------------------------
            if plotCoolantGraphFlag
                figure('Name', CoolantGroupName, 'Color','w');
                gp = plot(DirectedGraphCoolant);
                gp.Interpreter = 'none';
                title(sprintf('Directed Flow Graph: %s', CoolantGroupName), 'Interpreter','none');
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

        function [AdjacencyMat, InletArrayIdxsCell, OutletArrayIdxsCell, CoolantArrayIdxs] = getAdjacencyMatAndInletOutletIdxs(obj)
            % Return adjacency matrix, plus cell arrays of inlets/outlets for each enabled cooling system.
        
            IdxsEnabledCoolSys = find(obj.EnabledCoolingSystems);
            nEnabled = length(IdxsEnabledCoolSys);
        
            InletArrayIdxsCell = cell(nEnabled,1);
            OutletArrayIdxsCell = cell(nEnabled,1);
        
            CoolantArrayIdxs = [];
            AdjacencyMat = zeros(size(obj.ResMat));  % global adjacency
        
            for iSys = 1:nEnabled
                idxCoolSys = IdxsEnabledCoolSys(iSys);
                if strcmp(obj.SupportedCoolingSystemsMcadNames{idxCoolSys}, 'Blown Over')
                    % Blown Over is a special case: inlet and outlet is the Ambient node
                    ambientMcadIdx = 0; % in MotorCAD, ambient node always has index = 0
                    ambientArrayIdx = find([obj.NodeNamesAndMcadIdx{:,2}]==ambientMcadIdx);
                    InletArrayIdxsCell{iSys} = ambientArrayIdx;
                    OutletArrayIdxsCell{iSys} = ambientArrayIdx;
                end
                thisDigraphStruct = obj.CoolingSystemsDigraphs(idxCoolSys);
                thisDigraph = thisDigraphStruct.Digraph;
                if isempty(thisDigraph)
                    % system not enabled or no nodes
                    continue;
                end
        
                thisNodeNamesAndArrayIdx = thisDigraphStruct.NodeNamesAndArrayIdx;
                nodeNames      = thisNodeNamesAndArrayIdx(:,1);
                nodeArrayIdxs  = [thisNodeNamesAndArrayIdx{:,2}];
                localAdjMat    = adjacency(thisDigraph);
        
                % Add to the global adjacency
                AdjacencyMat(nodeArrayIdxs, nodeArrayIdxs) = localAdjMat;
        
                % Gather all these coolant indices
                CoolantArrayIdxs = [CoolantArrayIdxs; nodeArrayIdxs(:)]; %#ok<AGROW> 
        
                % Find inlets by name (matching how we do it above)
                isInlet = contains(nodeNames, 'Inlet', 'IgnoreCase', true);
                InletArrayIdxsCell{iSys} = nodeArrayIdxs(isInlet);
        
                % Find outlets by outdegree=0
                od = outdegree(thisDigraph);
                isOutlet = (od == 0) & ~isInlet; 
                OutletArrayIdxsCell{iSys} = nodeArrayIdxs(isOutlet);
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
                         coolingSystemsEnabled, BkptsStruct, options)
            % Generate a Simulink model implementing a ROM consisting of a 
            % set of state-space models at certain breakpoints (shaft speeds, 
            % coolant flow rates, and coolant inlet temperatures). The ROM
            % interpolates the state-space arrays between breakpoints.

            arguments
                obj
                modelName (1,1) string
                coolingSystemsEnabled (1,:) cell
                BkptsStruct (1,1) struct
                options.DCBusVoltage = []
            end
            
            assertBkptsStruct(BkptsStruct);

            obj.turnOffLossDependenceWithTemperatureOrSpeed();

            % 1) COMPUTE ROM DATA  ---------------------------------------
            
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
                if strcmp(thisCoolingSys, 'Blown Over')
                    frPropNames{idxCoolingSys} = 'BlownOver_FlowVelocity_mps';
                else
                    frPropNames{idxCoolingSys} = strcat(erase(thisCoolingSys, ' '), '_FlowRate_m3ps');
                end
                if any(startsWith(obj.MultiNozzleSprayCoolingSystems, string(thisCoolingSys))) % multi-nozzle cooling (Front and Rear)
                    obj.SprayCooling_Enable          = int32(1); % turn on the spray cooling system
                    obj.SprayCoolingNozzleDefinition = int32(1); % turn on multi-nozzle
                    TinPropNames{idxCoolingSys} = {...;
                        strcat(erase(thisCoolingSys, ' '), '_InletTemperature_F_degC');
                        strcat(erase(thisCoolingSys, ' '), '_InletTemperature_R_degC');
                        };
                else % normal case
                    TinPropNames{idxCoolingSys} = strcat(erase(thisCoolingSys, ' '), '_InletTemperature_degC');
                end
            end
            if numel(options.DCBusVoltage)==1
                obj.DCBusVoltage=options.DCBusVoltage;
            end
            obj.updateModel();
            AdjacencyMat = obj.AdjacencyMat;

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
            
            % Build inlet/outlet index rows that match coolingSystemsEnabled one-to-one
            numCoolSys   = numel(coolingSystemsEnabled);
            inletCells   = cell(1,numCoolSys);
            outletCells  = cell(1,numCoolSys);           
            for kSys = 1:numCoolSys
                csName = coolingSystemsEnabled{kSys};
                % identify which rows in the object belong to csName
                if any(startsWith(obj.MultiNozzleSprayCoolingSystems, string(csName))) % is multi-nozzle spray cooling (front and rear)
                    % Multi-nozzle aggregate: collect all spray rows
                    rowMask = ismember(obj.CoolingSystemNamesAndMcadIdxes(:,1), ...
                                       obj.MultiNozzleSprayCoolingSystems);
                else
                    % Normal case: exact match
                    rowMask = strcmp(obj.CoolingSystemNamesAndMcadIdxes(:,1), csName);
                end
                rows = find(rowMask);
                if isempty(rows)
                    error(['Cooling system "%s" is not enabled in Motor-CAD (or its ' ...
                           'name does not match SupportedCoolingSystemsMcadNames).'], csName);
                end
                % concatenate indices from all matching rows
                inIdx  = [];
                outIdx = [];
                for r = rows(:)'
                    inIdx  = [inIdx  obj.InletArrayIdxs{r}];           %#ok<AGROW>
                    tmpOut = obj.OutletArrayIdxs{r};
                    if isempty(tmpOut)          % per-element fallback
                        tmpOut = obj.InletArrayIdxs{r};
                    end
                    outIdx = [outIdx tmpOut];   %#ok<AGROW>
                end
                inletCells{kSys}  = inIdx;
                outletCells{kSys} = outIdx;
            end
            obj.InletArrayIdxs = inletCells;
            obj.OutletArrayIdxs = outletCells;

            % Simulink does not accept cell arrays as parameters, need to
            % transform them to arrays:
            InletCoolIdxs = adaptInletOutletArrays(obj.InletArrayIdxs); 
            OutletCoolIdxs = adaptInletOutletArrays(obj.OutletArrayIdxs);
            
            % Calculate state-space model at each breakpoint -------
            idxsCell = getNestedForLoopIdxs(sizeMatND(1:end-2));

            % Calculate state-space matrices for each breakpoint
            disp("Calculating state-space matrices for each breakpoint...");
            for idxBkptCmb = 1:length(idxsCell)
                bkptIdxsComb = idxsCell{idxBkptCmb};
                idxSpeed = bkptIdxsComb(1);
                disp("Breakpoint #" + num2str(idxBkptCmb) + " of " + num2str(length(idxsCell)));
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
                    thisCoolingSystem = coolingSystemsEnabled{idxCool};
                    if strcmp(thisCoolingSystem, 'Blown Over')
                        disp("  " + thisCoolingSystem + " flow velocity = " + frVal + " m/s");
                        obj.(frPropNames{idxCool}) = frVal;
                    else
                        disp("  " + thisCoolingSystem + " flow rate = " + frVal + " lpm");
                        obj.(frPropNames{idxCool}) = frVal/60/1000;
                    end
                    TinVal = TinBkpts(bkptIdxsComb(idxTin));
                    if iscell(TinPropNames{idxCool}) % multi-nozzle (F and R)
                        TinMultiNozzlePropNames = TinPropNames{idxCool};
                        disp("  " + thisCoolingSystem +" F inlet temperature = " + TinVal + " degC");
                        obj.(TinMultiNozzlePropNames{1}) = TinVal;
                        disp("  " + thisCoolingSystem +" R inlet temperature = " + TinVal + " degC");
                        obj.(TinMultiNozzlePropNames{2}) = TinVal;
                    else
                        disp("  " + thisCoolingSystem +" inlet temperature = " + TinVal + " degC");
                        obj.(TinPropNames{idxCool}) = TinVal;
                    end
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

            if numel(options.DCBusVoltage) > 1 % 3-D loss maps (LossMapsWithDCV)
                disp("Calculating losses at each DC Bus Voltage breakpoint ...")
                % Pre-allocate array sizes
                [ShaftTorqueVec, SpeedVec, ~] = reInterpolateTable( ...
                        obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Stator_Copper_Loss_Mat);
                DCBusVoltageVec = options.DCBusVoltage; 
                nX = numel(ShaftTorqueVec);           % Torque axis length
                nY = numel(SpeedVec);                 % Speed  axis length
                nZ = numel(DCBusVoltageVec);     % Voltage sweep length           
                % Each loss map has size [speed  × torque × voltage]
                xStator_Copper_Loss_Mat          = zeros(nX, nY, nZ);
                xRotor_Cage_Loss_Mat             = zeros(nX, nY, nZ);
                xIron_Loss_Stator_Back_Iron_Mat  = zeros(nX, nY, nZ);
                xIron_Loss_Stator_Tooth_Mat      = zeros(nX, nY, nZ);
                xStray_Load_Loss_Mat             = zeros(nX, nY, nZ);
                xMagnet_Loss_Mat                 = zeros(nX, nY, nZ);
                xIron_Loss_Rotor_Pole_Mat        = zeros(nX, nY, nZ);
                xIron_Loss_Rotor_Back_Iron_Mat   = zeros(nX, nY, nZ);
                xIron_Loss_Rotor_Tooth_Mat       = zeros(nX, nY, nZ);
                xFriction_Loss_Mat               = zeros(nX, nY, nZ);
                xWindage_Loss_Mat                = zeros(nX, nY, nZ);
                xStator_Copper_Loss_AC_Mat       = zeros(nX, nY, nZ);
                xBanding_Loss_Mat                = zeros(nX, nY, nZ);
                xSleeve_Loss_Mat                 = zeros(nX, nY, nZ);
            
                % Sweep the DC-bus voltages, refresh Motor-CAD, and stack the resulting 2-D maps into the pre-allocated 3-D arrays.
                for idxDCV = 1:nZ

                    disp("  DCBusVoltage = " + num2str(DCBusVoltageVec(idxDCV)) + " V");
                    % -- Set the DCV value in Motor-CAD
                    obj.DCBusVoltage =  DCBusVoltageVec(idxDCV);
            
                    % -- Ask Motor-CAD to recalculate loss tables at this voltage
                    obj.calculateMagneticLab(); % force lab .mat file refresh
                    obj.updateLabLossTables();
            
                    % -- Re-interpolate each refreshed loss table onto the common grid
                    [~,~,xStator_Copper_Loss_Mat(:,:,idxDCV)]         = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Stator_Copper_Loss_Mat);
                    [~,~,xRotor_Cage_Loss_Mat(:,:,idxDCV)]            = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Rotor_Cage_Loss_Mat);
                    [~,~,xIron_Loss_Stator_Back_Iron_Mat(:,:,idxDCV)] = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Stator_Back_Iron_Mat);
                    [~,~,xIron_Loss_Stator_Tooth_Mat(:,:,idxDCV)]     = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Stator_Tooth_Mat);
                    [~,~,xStray_Load_Loss_Mat(:,:,idxDCV)]            = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Stray_Load_Loss_Mat);
                    [~,~,xMagnet_Loss_Mat(:,:,idxDCV)]                = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Magnet_Loss_Mat);
                    [~,~,xIron_Loss_Rotor_Pole_Mat(:,:,idxDCV)]       = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Rotor_Pole_Mat);
                    [~,~,xIron_Loss_Rotor_Back_Iron_Mat(:,:,idxDCV)]  = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Rotor_Back_Iron_Mat);
                    [~,~,xIron_Loss_Rotor_Tooth_Mat(:,:,idxDCV)]      = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Iron_Loss_Rotor_Tooth_Mat);
                    [~,~,xFriction_Loss_Mat(:,:,idxDCV)]              = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Friction_Loss_Mat);
                    [~,~,xWindage_Loss_Mat(:,:,idxDCV)]               = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Windage_Loss_Mat);
                    [~,~,xStator_Copper_Loss_AC_Mat(:,:,idxDCV)]      = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Stator_Copper_Loss_AC_Mat);
                    [~,~,xBanding_Loss_Mat(:,:,idxDCV)]               = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Banding_Loss_Mat);
                    [~,~,xSleeve_Loss_Mat(:,:,idxDCV)]                = reInterpolateTable(obj.Shaft_Torque_Mat, obj.Speed_Mat, obj.Sleeve_Loss_Mat);
                end
            else % 2D loss maps
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
            end


            % Get PowerLossDistributor params
            disp("Estimating power loss distribution for each loss type...")
            lossDistrForEachType = obj.getLossDistrForEachType();
            TrefStator = obj.TrefStator_degC;
            TrefRotor = obj.TrefRotor_degC;
            xStatorCopperTempCoefResistivity = obj.StatorTempCoeffRes;
            xRotorCopperTempCoefResistivity = obj.RotorTempCoeffRes;
            xCoolantArrayIdxs = obj.CoolantArrayIdxs;            
            xCapMat = obj.CapMat;

            % 2) GENERATE SIMULINK MODEL ---------------------------------
            
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
            hOutTnodes = add_block('simulink/Commonly Used Blocks/Out1', strcat(romSubsysPath, '/TempNodes'));
            hInTorque = add_block('simulink/Commonly Used Blocks/In1', strcat(romSubsysPath, '/ShaftTorque_Nm'));
            hInSpeed = add_block('simulink/Commonly Used Blocks/In1', strcat(romSubsysPath, '/ShaftSpeed_RPM')); 
            if numel(options.DCBusVoltage) > 1 % 3-D loss maps (LossMapsWithDCV)
                hInVoltage = add_block('simulink/Commonly Used Blocks/In1', strcat(romSubsysPath, '/DCBusVoltage')); 
                hLossMaps = add_block(strcat(obj.mcadROMLibName, '/LossMapsWithDCV'), ...
                    strcat(romSubsysPath, '/Loss Maps from Motor-CAD Lab'), ...
                    'DCBusVoltageVec', 'DCBusVoltageVec', ...
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
            else  % 2-D loss maps (LossMapsWithDCV)
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
            end
            
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
            if numel(options.DCBusVoltage) > 1 % Voltage included
                hInVoltagePort = get_param(hInVoltage, 'PortHandles');
                add_line(romSubsysPath, hInVoltagePort.Outport, hLossMapsPorts.Inport(3));
            end
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
            if numel(options.DCBusVoltage) > 1 % Voltage included
                assignin(mdlWks,'DCBusVoltageVec', DCBusVoltageVec);
            end
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
            if numel(options.DCBusVoltage) > 1 % Voltage included
                add_block('simulink/Sources/Constant', strcat(modelName, '/DCBusVoltage'), ...
                'Value', num2str(options.DCBusVoltage(end)));
                add_line(modelName, strcat('DCBusVoltage', '/1'), strcat(motName, '_ReducedOrderModel', '/3'));
                coolant2ROMInputPortNumberOffset = 3;
            else
                coolant2ROMInputPortNumberOffset = 2;
            end
            add_line(modelName, strcat('TorqueNm', '/1'), strcat(motName, '_ReducedOrderModel', '/1'))
            add_line(modelName, strcat('SpeedRPM', '/1'), strcat(motName, '_ReducedOrderModel', '/2'))
            for idxCool = 1:numCoolSys
                thisCoolingSys = coolingSystemsEnabled{idxCool};
                frBlockName = strcat(erase(thisCoolingSys, ' '), '_FlowRate_lpm');
                TinBlockName = strcat(erase(thisCoolingSys, ' '), '_InletTemp_degC');
                add_block('simulink/Sources/Constant', strcat(modelName, '/', frBlockName), ...
                'Value', '3');
                add_block('simulink/Sources/Constant', strcat(modelName, '/', TinBlockName), ...
                'Value', '25');
                add_block('simulink/Commonly Used Blocks/Mux', strcat(modelName, '/Mux', num2str(idxCool)));
                add_line(modelName, strcat(frBlockName, '/1'), strcat('Mux', num2str(idxCool), '/1'))
                add_line(modelName, strcat(TinBlockName, '/1'), strcat('Mux', num2str(idxCool), '/2'))
                add_line(modelName, strcat('Mux', num2str(idxCool), '/1'), strcat(motName, '_ReducedOrderModel', '/', num2str(coolant2ROMInputPortNumberOffset+idxCool)));               
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

        function csData = splitSprayCoolingGroups(obj)
            % This helper returns an updated cell array that replaces the single
            % "Spray Cooling" row with:
            %   - A single "Spray Cooling" row in single-nozzle mode, or
            %   - Up to six subsystem rows in multi-nozzle mode (one for each of the
            %     front and rear nozzles for the enabled spray cooling systems).
            
            % Start by copying all group data except the 'Spray Cooling' row.
            allGroups = obj.GroupNamesAndMcadIdxes;
            isSprayCoolingGroup = strcmp(allGroups(:,1), 'Spray Cooling');
            csData = allGroups(~isSprayCoolingGroup,:);
            
            % If there is no "Spray Cooling" group, exit early.
            if ~any(isSprayCoolingGroup)
                return;
            end
            
            % Extract the spray cooling MCAD indexes.
            scIdx = find(isSprayCoolingGroup, 1, 'first');
            sprayCoolingMcadIdxes = allGroups{scIdx,2};  % e.g. [192 193 196 197]
            
            if obj.SprayCoolingNozzleDefinition == 0
                %=== CASE 1: Single-nozzle mode
                % Keep the entire "Spray Cooling" group as a single system.
                csData(end+1,:) = {'Spray Cooling', sprayCoolingMcadIdxes};
                
            else
                %=== CASE 2: Multi-nozzle mode
                % Each enabled spray cooling system will be split into two independent
                % cooling subsystems (front and rear).
                
                % Retrieve the node names corresponding to the spray cooling MCAD indexes.
                allMcadIdxes = [obj.NodeNamesAndMcadIdx{:,2}];
                allNodeNames = obj.NodeNames;  % cell array of all node names, in order
                
                % Find the positions in allNodeNames corresponding to sprayCoolingMcadIdxes.
                [~, loc] = ismember(sprayCoolingMcadIdxes, allMcadIdxes);
                nodeNamesInGroup = allNodeNames(loc);  % node names in the "Spray Cooling" group
                
                % Define the base spray cooling system names and their enable flags.
                % (There is one flag per system; each will yield two rows: _F and _R.)
                baseSysNames = {'Spray_RadialHousing', 'Spray_RadialRotor', 'Spray_AxialEndcap'};
                baseSysEnables = [obj.Spray_RadialHousing_Enable, obj.Spray_RadialRotor_Enable, obj.Spray_AxialEndcap_Enable];
                
                % Define the two nozzle sides.
                sides = {'_F', '_R'};
                
                % Loop over each base system.
                for iSys = 1:length(baseSysNames)
                    if baseSysEnables(iSys)
                        for iSide = 1:length(sides)
                            % Construct the expected node name substring for this subsystem.
                            subSysName = [baseSysNames{iSys} sides{iSide}];
                            
                            % Identify the inlet index corresponding to this subSysName
                            isThisSubSys = startsWith(nodeNamesInGroup, baseSysNames{iSys}) ...
                                                & endsWith  (nodeNamesInGroup, sides{iSide});
                            theseNodeMcadIdx = sprayCoolingMcadIdxes(isThisSubSys);

                            % Add a new row for this cooling subsystem.
                            csData(end+1,:) = {subSysName, theseNodeMcadIdx}; %#ok<AGROW> 
                        end
                    end
                end
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

function [array1d, mcadIdxes] = readMfFile1D(FileName)
    %READMFFILE1D Reads a 1D Motor-CAD matrix file (.tmf, .cmf, .pmf) and
    %extracts values and associated Motor-CAD indices from it.
    %
    % [array1d, mcadIdxes] = readMfFile1D(FileName)
    %
    % Inputs:
    %   FileName - String. Path to the .tmf, .cmf, or .pmf file.
    %
    % Outputs:
    %   array1d   - Numeric array. Values associated with each node.
    %   mcadIdxes - Numeric array. Corresponding Motor-CAD indices of the nodes.

    fid = fopen(FileName, 'rt');
    if fid < 0
        error('Cannot open the file: %s', FileName);
    end

    % Skip header lines
    headerLine1 = fgetl(fid); %#ok<NASGU>
    headerLine2 = fgetl(fid); %#ok<NASGU>
    headerLine3 = fgetl(fid); %#ok<NASGU> % Possibly the "Number of nodes" line

    % Initialize containers
    array1d = [];
    mcadIdxes = [];

    % Updated regular expression to handle lines with or without closing parenthesis
    %
    % Regex Explanation:
    % -----------------------------------
    % ^\s*                           - Match any leading whitespace at the start of the line
    % (\d+)                          - Capture the node index (one or more digits) as Group 1
    % .*?                            - Lazily match any characters (node name), allowing for incomplete names
    % ([-]?\d+\.?\d*(?:[Ee][-+]?\d+)?) - Capture the numeric value (temperature) as Group 2
    % \s*;?\s*$                      - Match optional whitespace and an optional semicolon until the end of the line
    expr = '^\s*(\d+).*?([-]?\d+\.?\d*(?:[Ee][-+]?\d+)?)\s*;?\s*$';

    lineNumber = 4; % Starting from the 4th line after headers

    while true
        thisLine = fgetl(fid);
        if ~ischar(thisLine)
            % End of file
            break;
        end
        thisLine = strtrim(thisLine);
        if isempty(thisLine)
            % Skip empty lines
            continue;
        end
        tokens = regexp(thisLine, expr, 'tokens', 'once');
        if ~isempty(tokens)
            idxStr = tokens{1};
            valStr = tokens{2};

            % Convert strings to appropriate types
            mcadIdx = str2double(idxStr);
            val = str2double(valStr);

            % Append to arrays
            mcadIdxes(end+1) = mcadIdx; %#ok<AGROW>
            array1d(end+1) = val; %#ok<AGROW>
        else
            warning('Line %d in "%s" did not match expected pattern and was skipped.', lineNumber, FileName);
        end
        lineNumber = lineNumber + 1;
    end

    fclose(fid);
end

function array2d = readMfFile2D(filename, numNodes)
    fid = fopen(filename, 'r');
    if fid < 0
        error("Could not open file '%s'.", filename);
    end

    %-----------------------------------------------------------
    % 1) Skip any fixed header lines you do not need
    %-----------------------------------------------------------
    for i = 1:4
        fgetl(fid); % just discard these lines
    end

    % Initialize your matrix
    array2d = zeros(numNodes, numNodes);

    %-----------------------------------------------------------
    % 2) Read line by line for 'numNodes' lines
    %-----------------------------------------------------------
    for row = 1:numNodes
        line = fgetl(fid);
        if ~ischar(line)
            % If we hit end-of-file before reading numNodes lines
            error('Unexpected end of file at line %d.', row + 3);
        end

        %-------------------------------------------------------
        % 3) Split the line by semicolons
        %    (Adjust delimiter if needed; e.g. commas, spaces, etc.)
        %-------------------------------------------------------
        parts = strsplit(line, ';');
        
        % If the last split is empty (common when lines end in ";"), remove it.
        if isempty(parts{end})
            parts(end) = [];
        end

        %-------------------------------------------------------
        % 4) Convert each piece to a number when possible
        %-------------------------------------------------------
        numericVals = [];
        for k = 1:numel(parts)
            val = str2double(strtrim(parts{k}));
            if ~isnan(val)
                numericVals(end+1) = val; %#ok<AGROW>
            end
        end

        %-------------------------------------------------------
        % 5) Make sure we got exactly 'numNodes' numeric values
        %    for this line. Adjust as needed for your file format.
        %-------------------------------------------------------
        if length(numericVals) < numNodes
            % Probably the first chunk includes a truncated node name and
            % the first numeric value (e.g. "192 (Spray_RadialHousing_Inlet_Flui 1000000000...")
            % We need to extract that extra number and prepend it to numericVals.
        
            firstChunk = parts{1};  
            % Use a regex to find all numbers in that chunk:
            tokens = regexp(firstChunk,'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)','match'); % matches numbers with int, float, and scientific format
            % tokens might look like {"192","1000000000"}.
        
            if length(tokens) >= 2
                % The second number is the data (first number is the Motor-CAD index, so skip)
                firstElement = str2double(tokens{2});
                % Prepend it to our existing numeric values
                numericVals = [firstElement, numericVals]; %#ok<AGROW> 
            else
                fclose(fid);
                error(['Could not recover a missing numeric value from line %d. ' ...
                       'Check for truncation or missing delimiters.'], row + 3);
            end
        
            % After fixing, make sure we now have enough data:
            if length(numericVals) < numNodes
                fclose(fid);
                error('Still not enough data on line %d, even after attempted fix.', row + 3);
            end
        
        elseif length(numericVals) > numNodes
            fclose(fid);
            error('Too many numeric elements on line %d. Check file format.', row + 3);
        end

        % Assign row data to the output matrix
        array2d(row, :) = numericVals;
    end

    fclose(fid);
end

function [GroupNames, NodeNames, GroupIdxs, McadIdxs] = readNmfFile(FileName)
    %READNMFFILE Reads a Node Grouping Matrix (.nmf) file and extracts group names,
    %node names, and their corresponding Motor-CAD indices in a robust manner.
    %
    % [GroupNames, NodeNames, GroupIdxs, McadIdxs] = readNmfFile(FileName)
    %
    % Inputs:
    %   FileName - String. Path to the .nmf file.
    %
    % Outputs:
    %   GroupNames - Cell array of strings. Names of the groups.
    %   NodeNames  - Cell array of strings. Names of the nodes.
    %   GroupIdxs  - Cell array of numeric arrays. Motor-CAD indices of nodes in each group.
    %   McadIdxs   - Cell array of doubles. Motor-CAD index of each node.

    fid = fopen(FileName, 'rt');
    if fid < 0
        error('Cannot open the file: %s', FileName);
    end

    % Skip header lines
    headerLine1 = fgetl(fid); %#ok<NASGU> 
    headerLine2 = fgetl(fid); %#ok<NASGU> 
    headerLine3 = fgetl(fid); %#ok<NASGU> 

    GroupNames = {};
    GroupIdxs  = {};
    NodeNames  = {};
    McadIdxs    = {};

    currentGroup = '';
    currentGroupIdxs = [];

    while true
        thisLine = fgetl(fid);
        if ~ischar(thisLine)
            % End of file
            break;
        end
        thisLine = strtrim(thisLine);
        if isempty(thisLine)
            % Skip empty lines
            continue;
        end
        if startsWith(thisLine, '[') && endsWith(thisLine, ']')
            % New group definition, e.g., [Armature Winding (Active)]
            % Store the previous group if it exists
            if ~isempty(currentGroup)
                GroupNames{end+1} = currentGroup; %#ok<AGROW>
                GroupIdxs{end+1}  = currentGroupIdxs; %#ok<AGROW>
            end
            % Start a new group
            currentGroup = thisLine(2:end-1);  % Remove the [ and ]
            currentGroupIdxs = [];
        else
            % Expected format: "Index (Name)" e.g., "342 (Wedge)"
            expr = '^(\d+)\s*\((.*)\)$';
            tokens = regexp(thisLine, expr, 'tokens', 'once');
            if ~isempty(tokens)
                idxStr    = tokens{1};
                nameStr   = tokens{2};
                idxVal    = str2double(idxStr);
                nodeName  = stripOuterParentheses(nameStr);

                NodeNames{end+1} = nodeName; %#ok<AGROW>
                McadIdxs{end+1}   = idxVal;   %#ok<AGROW>
                currentGroupIdxs(end+1) = idxVal; %#ok<AGROW>
            else
                warning('Line "%s" did not match expected pattern in %s', thisLine, FileName);
            end
        end
    end

    % Store the last group
    if ~isempty(currentGroup)
        GroupNames{end+1} = currentGroup;
        GroupIdxs{end+1}  = currentGroupIdxs;
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
        % Re-interpolate 2D scattered map into square grid

        warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId'); % suppress this unimportant warning temporarily
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
        warning('on', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId'); % re-enable the warning

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

function nodeNameOut = stripOuterParentheses(nodeNameIn)
    %STRIPOUTERPARENTHESES Removes the outermost parentheses from a node name if present.
    %
    % nodeNameOut = stripOuterParentheses(nodeNameIn) returns the node name
    % with the outermost '(' and ')' removed if they exist. Otherwise, returns
    % the original node name unchanged.
    
    nodeNameIn = strtrim(nodeNameIn);  % Remove leading/trailing whitespace
    if startsWith(nodeNameIn, '(') && endsWith(nodeNameIn, ')')
        % Remove the very first '(' and the very last ')'
        nodeNameOut = nodeNameIn(2:end-1);
    else
        nodeNameOut = nodeNameIn;
    end
end

function a = adaptInletOutletArrays(c)
    arguments
        c (1,:) cell
    end
    %   A = adaptInletOutletArrays(C) takes a 1×N cell array C whose elements
    %   are numeric vectors (row- or column-oriented) and returns a MxN
    %   double array where M is the length of the longest element in C and
    %   the rows that don't have enough elements have repeated the first
    %   element of the row to fit the MxN size
    %
    %   Example:
    %     c = {[0]; [1,2,3]; [4;5;6]; [7 8]};
    %     a = adaptInletOutletArrays(c);
    %     % a is [0 0 0; 1 2 3; 4 5 6; 7 8 7]

    % Find the length of each vector and the overall maximum length (M)
    vecLen = cellfun(@numel, c);
    M      = max(vecLen);        % longest vector length
    N      = numel(c);           % number of vectors / rows in the output

    % Pre-allocate the result (each input vector becomes one row)
    a = zeros(N, M);

    % Build each row, padding with the vector’s first element when needed
    for k = 1:N
        v = c{k}(:).';           % force row orientation
        if numel(v) < M
            v = [v repmat(v(1), 1, M - numel(v))]; %#ok<AGROW> 
        end
        a(k, :) = v;
    end
end

function assertBkptsStruct(BkptsStruct)
%ASSERTBKPTSSTRUCT  Validate breakpoint-definition struct.
%
% Rules enforced
% ---------------
% 1) BkptsStruct must be a scalar struct.
% 2) Required fields:  w, fr1, Tin1
% 3) Any other field must match  frN  or  TinN  where N is an integer ≥ 2.
% 4) Every frN must have a matching TinN with the same N (and vice-versa).
% 5) Every field’s value must contain at least two elements.
%
% Throws a descriptive error if any check fails.

    %-- Basic struct check
    if ~isstruct(BkptsStruct) || numel(BkptsStruct) ~= 1
        error('BkptsStruct must be a scalar struct.');
    end

    fn = fieldnames(BkptsStruct);

    %-- (2) Required fields present?
    req = {'w','fr1','Tin1'};
    missing = setdiff(req, fn);
    if ~isempty(missing)
        error('BkptsStruct is missing required field(s): %s', strjoin(missing, ', '));
    end

    %-- (3) / (4) Validate names and pairing
    frNums  = [];   % numeric suffixes for frN fields
    tinNums = [];   % numeric suffixes for TinN fields
    for k = 1:numel(fn)
        name = fn{k};

        % Allow only: 'w', 'frN', 'TinN'
        if strcmp(name,'w')
            continue
        end

        m = regexp(name,'^(fr|Tin)(\d+)$','tokens','once');
        if isempty(m)
            error('Invalid field name "%s" for BkptsStruct. Allowed names: w, frN, TinN', name);
        end

        suffix = str2double(m{2});
        if suffix < 1 || floor(suffix) ~= suffix
            error('BkptsStruct field "%s" has invalid suffix (must be an integer ≥ 1).', name);
        end

        if strcmp(m{1},'fr')
            frNums(end+1)  = suffix;  %#ok<AGROW>
        else
            tinNums(end+1) = suffix;  %#ok<AGROW>
        end
    end

    % Make sure frN/TinN pairs match exactly
    if ~isequal(sort(frNums), sort(tinNums))
        error('Invalid BkptsStruct. Each "frN" field must have a matching "TinN" field (and vice-versa).');
    end

    %-- (1) / (5) Check that each field has ≥ 2 values
    for k = 1:numel(fn)
        v = BkptsStruct.(fn{k});
        if numel(v) < 2
            error('Invalid BkptsStruct. Field "%s" must contain at least two values.', fn{k});
        end
    end
end
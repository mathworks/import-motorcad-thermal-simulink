classdef (Abstract) AbstractInterface < handle
    %ABSTRACTINTERFACE Base class for Motor-CAD interfaces that implements
    % the basic load, set, and get methods.
    % The client interface classes require a mapping of Motor-CAD parameter
    % names to MATLAB parameter names.
    
    % Copyright 2022-2023 The MathWorks, Inc.
    
    properties(Abstract, Constant, Access=protected)
        McadParameterNameList cell % Motor-CAD parameter names list
        ParameterNameList cell % Interface class parameter names list
    end

    properties(Access=protected)
        mcad % Motor-CAD instance object
        motFullFile % mot file full path
        McadParameterValuesList cell % Motor-CAD parameter values list
    end
    
    methods(Access=public)
        function obj = AbstractInterface(motFullFile)
            % AbstractInterface constructor method
            obj.loadMcad(motFullFile);
            obj.getMcadParameters();
            obj.setMcadMatrixSeparator();

            for idxParam = 1:length(obj.McadParameterNameList)
                paramName = obj.ParameterNameList{idxParam};
                paramValue = obj.McadParameterValuesList{idxParam};
                obj.(paramName) = paramValue;
            end
        end

        function setParameter(obj, paramName, paramValue)
            % Set parameter value in Motor-CAD
            idxParam = strcmp(obj.ParameterNameList, paramName);
            obj.McadParameterValuesList{idxParam} = paramValue;
            mcadParameterName = obj.McadParameterNameList{idxParam};
            invoke(obj.mcad,'DisplayScreen','Scripting');
            invoke(obj.mcad, 'SetVariable', mcadParameterName, paramValue);
        end
        
        function paramValue = getParameter(obj,paramName)
            % Get parameter value from Motor-CAD
            idxParam = strcmp(obj.ParameterNameList, paramName);
            paramValue = obj.McadParameterValuesList{idxParam};
        end
    end

    methods(Access=private)
        function loadMcad(obj, motFullFile)  
            % Load the Motor-CAD .mot file
            motFullFile = which(motFullFile);
            obj.motFullFile = motFullFile;
            if ~isfile(motFullFile)
                error('.mot file not found on the path. Add the file to the path or use the full file path.');
            end
            try
                obj.mcad = actxserver('MotorCAD.AppAutomation');
            catch theException
                % Construct a combined error message
                combinedErrorMessage = "MATLAB is unable to connect to the Motor-CAD ActiveX automation server. " + ...
                                               "Ensure you opened the .mot file with 'Visual Basic' selected " + ...
                                               "as Scripting Engine option." + newline + ...
                                               "Root cause error: " + newline + ...
                                               "Error ID: " + theException.identifier + newline + ...
                                               "Error Message: " + theException.message;
                % Throw a new error with the combined message
                error(combinedErrorMessage);
            end
            invoke(obj.mcad, 'SetVariable',"MessageDisplayState",2); % avoid pop-up windows
            invoke(obj.mcad,'LoadFromFile',motFullFile);
            invoke(obj.mcad,'DisplayScreen','Scripting');
        end

        function getMcadParameters(obj)
            % Get all Motor-CAD parameters in the McadParameterNameList
            for idxParam = 1:length(obj.McadParameterNameList)
                paramName = obj.McadParameterNameList{idxParam};
                [success, paramValue] = invoke(obj.mcad, 'GetVariable', paramName);
                if success ~= 0
                    warn('GetVariable failed')
                    paramValue = NaN;
                end
                obj.McadParameterValuesList{idxParam} = paramValue;
            end
        end

        function setMcadMatrixSeparator(obj)
            % Set the Motor-CAD MatrixTextSeparator property to ";"
            invoke(obj.mcad, 'SetVariable',"MatrixTextSeparator",";");
        end
    end
end
classdef (Abstract) AbstractInterface < handle
    %ABSTRACTINTERFACE Base class for Motor-CAD interfaces that implements
    % the basic load, set, and get methods.
    % The client interface classes require a mapping of Motor-CAD parameter
    % names to MATLAB parameter names.
    
    % Copyright 2022-2024 The MathWorks, Inc.
    
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
            obj.mcad.display_screen('Scripting');
            obj.mcad.set_variable(mcadParameterName, paramValue);
        end
        
        function paramValue = getParameter(obj,paramName)
            % Get parameter value from Motor-CAD
            idxParam = strcmp(obj.ParameterNameList, paramName);
            paramValue = obj.McadParameterValuesList{idxParam};
        end

        function delete(obj)
            % Destructor method to close the Motor-CAD instance before
            % deleting the object
            if ~isempty(obj.mcad)
                try
                    obj.mcad.quit();
                catch
                    warning('Failed to close the Motor-CAD instance on cleanup.');
                end
            end
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
                pymotorcad = py.importlib.import_module('ansys.motorcad.core');
            catch theException
                % Define a clickable hyperlink for the "PyMotorCAD scripting in MATLAB" page
                url1 = "https://motorcad.docs.pyansys.com/version/stable/user_guide/matlab_scripting.html";
                linkText1 = "PyMotorCAD scripting in MATLAB";
                % Define a clickable hyperlink for the "Versions of Python Compatible with MATLAB Products by Release" page
                url2 = "https://www.mathworks.com/support/requirements/python-compatibility.html";
                linkText2 = "Versions of Python Compatible with MATLAB Products by Release";
            
                % Construct the combined error message with a clickable hyperlink
                combinedErrorMessage = ...
                    "MATLAB is unable to connect to the PyMotorCAD interface. " + ...
                    "Ensure a compatible Python version is installed and " + ...
                    "the ansys.motorcad.core module is installed." + newline + ...
                    "For more information, see " + ...
                    "<a href=""" + url1 + """>" + linkText1 + "</a>" + " and " + ...
                    "<a href=""" + url2 + """>" + linkText2 + "</a>." + newline + ...
                    "Root cause error:" + newline + ...
                    "Error ID: " + theException.identifier + newline + ...
                    "Error Message: " + theException.message;
            
                % Throw a new error with the combined message
                error(combinedErrorMessage);
            end

            obj.mcad = pymotorcad.MotorCAD();
            obj.mcad.set_variable('MessageDisplayState',2); % avoid pop-up windows
            obj.mcad.load_from_file(motFullFile);
            obj.mcad.display_screen('Scripting');
        end

        function getMcadParameters(obj)
            % Get all Motor-CAD parameters in the McadParameterNameList
            for idxParam = 1:length(obj.McadParameterNameList)
                paramName = obj.McadParameterNameList{idxParam};
                paramValue = obj.mcad.get_variable(paramName);
                obj.McadParameterValuesList{idxParam} = double(paramValue);
            end
        end

        function setMcadMatrixSeparator(obj)
            % Set the Motor-CAD MatrixTextSeparator property to ";"
            obj.mcad.set_variable("MatrixTextSeparator",";"); % note: named ExportTextSeparator in v2024.2
        end
    end
end
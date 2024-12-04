function modelNamePath = findModelInProject(modelName)
    % FINDMODELINPROJECT returns the file path of a Simulink® model name if it
    % finds it in the current project up to two levels of directory depth
    % from the root project folder. Else, it returns an empty modelNamePath.
    % Input arguments:
    % - modelName: [string/char]: model name to search
    % Output arguments:
    % - modelNamePath: [char]: file path of the model name, if found, or
    % empty array, if not found.

    % Copyright 2022 The MathWorks, Inc.

    modelNamePath = [];
    theProject = matlab.project.rootProject;
    projectRootDir = theProject.RootFolder;
    dirStructRoot = dir(fullfile(projectRootDir, '*.slx'));
    dirStructOneLevelDeep = dir(fullfile(projectRootDir, '*', '*.slx'));
    dirStruct = [dirStructRoot; dirStructOneLevelDeep];  
    for idxFile = 1:length(dirStruct)
        fileName = dirStruct(idxFile).name;
        if strcmp(fileName, strcat(modelName, '.slx'))
            modelNamePath = fullfile(dirStruct(idxFile).folder, fileName);
            return
        end
    end
end


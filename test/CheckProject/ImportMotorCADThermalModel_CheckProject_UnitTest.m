classdef ImportMotorCADThermalModel_CheckProject_UnitTest < matlab.unittest.TestCase
% Class implementation of unit test

% Copyright 2023 The MathWorks, Inc.

methods (Test)

function CheckProject(~)
  close all
  bdclose all
  ImportMotorCADThermalModel_CheckProject
  close all
  bdclose all
end

end  % methods (Test)
end  % classdef

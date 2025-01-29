%% Script to run library unit tests
% This script runs all the unit tests that are the child classes of
% matlab.unittest.TestCase in the project.
% Unit test classes are automatically detected by
% the matlab.unittest.TestSuite.fromFolder function.

% Copyright 2022-2025 The MathWorks, Inc.

relstr = matlabRelease().Release;
disp("This is MATLAB " + relstr)

%% Create test suite

prjroot = currentProject().RootFolder;

suite = matlab.unittest.TestSuite.fromFolder(fullfile(prjroot, "test", "models"), "IncludingSubfolders",true);

%% Create test runner

runner = matlab.unittest.TestRunner.withTextOutput( ...
          "OutputDetail", matlab.unittest.Verbosity.Detailed);

%% JUnit style test result

plugin = matlab.unittest.plugins.XMLPlugin.producingJUnitFormat( ...
          fullfile(prjroot, "test", "ExampleModelTestResults_"+relstr+".xml"));

addPlugin(runner, plugin)

%% Run tests

results = run(runner, suite);
assertSuccess(results);
disp(results)

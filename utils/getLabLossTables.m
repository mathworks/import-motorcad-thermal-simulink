function [Speed, Shaft_Torque, ...
    Stator_Copper_Loss, ...
    Iron_Loss_Stator_Back_Iron, ...
    Iron_Loss_Stator_Tooth, ...
    Magnet_Loss, ...
    Iron_Loss_Rotor_Pole, ...
    Iron_Loss_Rotor_Back_Iron, ...
    Friction_Loss, ...
    Windage_Loss] = getLabLossTables(motFullFileName)

    % Copyright 2022 The MathWorks, Inc.
    
    [~, motFile] = fileparts(motFullFileName);
    LabFileName = fullfile(extractBefore(motFullFileName, motFile), motFile, 'Lab', 'MotorLAB_elecdata.mat');
    outLoad = load(LabFileName);

    Speed = outLoad.Speed;
    Shaft_Torque = outLoad.Shaft_Torque;
    Stator_Copper_Loss = outLoad.Stator_Copper_Loss;
    Iron_Loss_Stator_Back_Iron = outLoad.Iron_Loss_Stator_Back_Iron;
    Iron_Loss_Stator_Tooth = outLoad.Iron_Loss_Stator_Tooth;
    Magnet_Loss = outLoad.Magnet_Loss;
    Iron_Loss_Rotor_Pole = outLoad.Iron_Loss_Rotor_Pole;
    Iron_Loss_Rotor_Back_Iron = outLoad.Iron_Loss_Rotor_Back_Iron;
    Friction_Loss = outLoad.Friction_Loss;
    Windage_Loss = outLoad.Windage_Loss;

   
end

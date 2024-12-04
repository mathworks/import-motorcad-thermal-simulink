classdef BasicInterface < mcadinterface.AbstractInterface
    %BASICINTERFACE Motor-CAD interface that sets and gets parameter values and
    %implements callbacks to push changes of MATLAB parameter values to the
    %associated Motor-CAD instance.

    % Copyright 2022 The MathWorks, Inc.
    
    properties(Constant, Access=protected)
        McadParameterNameList = {'Ambient_Temperature'; ...
                                 'StatorCopperLossesVaryWithTemp'; ...
                                 'StatorCopperTempCoefResistivity'; ...
                                 'StatorWindingTemperatureAtWhichPcuInput'; ...
                                 'RotorCopperLossesVaryWithTemp'; ...
                                 'RotorCopperTempCoefResistivity'; ...
                                 'RotorWindingTemperatureAtWhichPcuInput'; ...
                                 'MagneticThermalCoupling'; ...
                                 'LabThermalCoupling'; ...
                                 ... Loss variables
                                 'Armature_Copper_Loss_@Ref_Speed'; ...
                                 'Rotor_Copper_Loss_@Ref_Speed'; ...
                                 'Stator_Iron_Loss_@Ref_Speed_[Back_Iron]'; ...
                                 'Stator_Iron_Loss_@Ref_Speed_[Tooth]'; ...
                                 'IM_SllStatorSplit_MotorLAB'; ...
                                 'Stator_Iron_Stray_Load_Loss_@Ref_Speed'; ...
                                 'Rotor_Iron_Stray_Load_Loss_@Ref_Speed'; ...
                                 'Stator_Copper_Stray_Load_Loss_@Ref_Speed'; ...
                                 'Rotor_Copper_Stray_Load_Loss_@Ref_Speed'; ...
                                 'Magnet_Iron_Loss_@Ref_Speed'; ...
                                 'Rotor_Iron_Loss_@Ref_Speed_[Embedded_Magnet_Pole]'; ...
                                 'Rotor_Iron_Loss_@Ref_Speed_[Back_Iron]'; ...
                                 'Rotor_Iron_Loss_@Ref_Speed_[Tooth]'; ...
                                 'Friction_Loss_[F]_@Ref_Speed'; ...
                                 'Friction_Loss_[R]_@Ref_Speed'; ...
                                 'Windage_Loss_@Ref_Speed'; ...
                                 'Windage_Loss_(Ext_Fan)@Ref_Speed'; ...
                                 ... 
                                 'TorqueDemand_MotorLAB'; ...
                                 'SpeedDemand_MotorLAB'; ...
                                 'Shaft_Speed_[RPM]'; ... 
                                 ... Additional loss options
                                 'Loss_Function_Speed'; ...
                                 'Copper_Losses_Vary_With_Temperature'; ...
                                 'RotorCopperLossesVaryWithTemp'; ...
                                 'StatorIronStrayLoadLossesVaryWithTemp'; ...
                                 'RotorIronStrayLoadLossesVaryWithTemp'; ...
                                 'StatorCopperStrayLoadLossesVaryWithTemp'; ...
                                 'RotorCopperStrayLoadLossesVaryWithTemp'; ...
                                 ... Cooling system variables -------
                                 ... Ventilated
                                 'Through_Ventilation'; ...
                                 'SelfVentilation'; ...
                                 'TVent_Flow_Rate'; ...
                                 'TVent_Inlet_Temperature'; ...
                                 ... Housing Water Jacket
                                 'Housing_Water_Jacket'; ...
                                 'WJ_Fluid_Volume_Flow_Rate'; ...
                                 'WJ_Fluid_Inlet_Temperature'; ...
                                 ... Shaft Spiral Groove
                                 'Shaft_Spiral_Groove'; ...
                                 'Shaft_Groove_Fluid_Volume_Flow_Rate'; ...
                                 'Shaft_Groove_Fluid_Inlet_Temperature'; ...
                                 ... Wet Rotor
                                 'Wet_Rotor'; ...
                                 'Wet_Rotor_Fluid_Volume_Flow_Rate'; ...
                                 'Wet_Rotor_Inlet_Temperature'; ...
                                 ... Spray Cooling
                                 'Spray_Cooling'; ...
                                 'Spray_Cooling_Fluid_Volume_Flow_Rate'; ...
                                 'Spray_Cooling_Inlet_Temp'; ...
                                 ... Rotor Water Jacket
                                 'Rotor_Water_Jacket'; ...
                                 'Rotor_WJ_Fluid_Volume_Flow_Rate'; ...
                                 'RotorWJ_Inlet_Temp'; ...
                                 ... Slot Water Jacket
                                 'Slot_Water_Jacket'; ...
                                 'Slot_WJ_Fluid_Volume_Flow_Rate'; ...
                                 'Slot_WJ_fluid_inlet_temperature'; ...
                                 ... add all cooling inlet variables
                                 ... Calculation settings
                                 'ThermalCalcType'; ...
                                 'TransientCalculationType'; ...
                                 'Simple_Transient_Period'; ...
                                 'EndPoint_Definition'; ...
                                 'Simple_Transient_Definition'; ...
                                 'Simple_Transient_Number_Points'; ...
                                 'Simple_Transient_Torque'; ...
                                 'InitialTransientTemperatureOption' ...
                                 };
        ParameterNameList = {'Tambient_degC'; ...
                             'EnableStatorTempCoeffRes'; ...
                             'StatorTempCoeffRes'; ...
                             'TrefStator_degC'; ...
                             'EnableRotorTempCoeffRes'; ...
                             'RotorTempCoeffRes'; ...
                             'TrefRotor_degC'; ...
                             'MagneticThermalCouplingChoice'; ...
                             'LabThermalCouplingChoice'; ...
                              ... Loss variables
                             'Armature_Copper_Loss'; ...
                             'Rotor_Copper_Loss'; ...
                             'Stator_Iron_Loss_Back_Iron'; ...
                             'Stator_Iron_Loss_Tooth'; ...
                             'Stray_Loss_Stator_Iron_Proportion'; ...
                             'Stator_Iron_Stray_Load_Loss'; ...
                             'Rotor_Iron_Stray_Load_Loss'; ...
                             'Stator_Copper_Stray_Load_Loss'; ...
                             'Rotor_Copper_Stray_Load_Loss'; ...
                             'Magnet_Iron_Loss'; ...
                             'Rotor_Iron_Loss_Embedded_Magnet_Pole'; ...
                             'Rotor_Iron_Loss_Back_Iron'; ...
                             'Rotor_Iron_Loss_Tooth'; ...
                             'Friction_Loss_F'; ...
                             'Friction_Loss_R'; ...
                             'Windage_Loss'; ...
                             'Windage_Loss_Ext_Fan'; ...
                             ...
                             'Shaft_Torque_Nm'; ...
                             'Shaft_Speed_RPM'; ...
                             'Shaft_Speed_RPM_Thermal'; ...
                             ... Additional loss options
                             'Loss_Function_Speed'; ...
                             'Copper_Losses_Vary_With_Temperature'; ...
                             'RotorCopperLossesVaryWithTemp'; ...
                             'StatorIronStrayLoadLossesVaryWithTemp'; ...
                             'RotorIronStrayLoadLossesVaryWithTemp'; ...
                             'StatorCopperStrayLoadLossesVaryWithTemp'; ...
                             'RotorCopperStrayLoadLossesVaryWithTemp'; ...
                             ... Cooling system variables -----------
                             ... Ventilated
                             'Ventilated_Enable'; ...
                             'Ventilated_Enable2'; ...
                             'Ventilated_FlowRate_m3ps'; ...
                             'Ventilated_InletTemperature_degC'; ...
                             ... Housing Water Jacket
                             'HousingWaterJacket_Enable'; ...
                             'HousingWaterJacket_FlowRate_m3ps'; ...
                             'HousingWaterJacket_InletTemperature_degC'; ...
                             ... Shaft Spiral Groove
                             'ShaftSpiralGroove_Enable'; ...
                             'ShaftSpiralGroove_FlowRate_m3ps'; ...
                             'ShaftSpiralGroove_InletTemperature_degC'; ...
                             ... Wet Rotor
                             'WetRotor_Enable'; ...
                             'WetRotor_FlowRate_m3ps'; ...
                             'WetRotor_InletTemperature_degC'; ...
                             ... Spray Cooling
                             'SprayCooling_Enable'; ...
                             'SprayCooling_FlowRate_m3ps'; ...
                             'SprayCooling_InletTemperature_degC'; ...
                             ... Rotor Water Jacket
                             'RotorWaterJacket_Enable'; ...
                             'RotorWaterJacket_FlowRate_m3ps'; ...
                             'RotorWaterJacket_InletTemperature_degC'; ...
                             ... Slot Water Jacket
                             'SlotWaterJacket_Enable'; ...
                             'SlotWaterJacket_FlowRate_m3ps'; ...
                             'SlotWaterJacket_InletTemperature_degC'; ...
                             ... add all cooling inlet variables
                             ... Calculation settings
                             'ThermalSteadyOrTransientChoice'; ...
                             'TransientOption'; ...
                             'StopTime'; ...
                             'StopTimeChoice'; ...
                             'TransientLossesOrTorqueChoice'; ...
                             'NumTimeSteps'; ...
                             'TransientTorqueValue'; ...
                             'InitialTemperatureOption' ...
                             };
    end

    properties(Access=public)
        Tambient_degC (1,1) double % Ambient temperature [degC]
        EnableStatorTempCoeffRes (1,1) int32 % Enable stator temperature coefficient for resistivity
        StatorTempCoeffRes (1,1) double % Stator temperature coefficient for resistivity [S.I.]
        TrefStator_degC (1,1) double % Stator reference temperature [degC]
        EnableRotorTempCoeffRes (1,1) int32 % Enable rotor temperature coefficient for resistivity
        RotorTempCoeffRes (1,1) double % Rotor temperature coefficient for resistivity [S.I.]
        TrefRotor_degC (1,1) double % Rotor reference temperature [degC]
        MagneticThermalCouplingChoice (1,1) int32 % Magnetic-Thermal coupling choice
        LabThermalCouplingChoice (1,1) int32 % Lab-Thermal coupling choice
        ...
        Shaft_Torque_Nm (1,1) double % Shaft torque operating point [Nm]
        Shaft_Speed_RPM (1,1) double % Shaft speed operating point [rpm]
        ... Additional loss options
        Loss_Function_Speed (1,1) int32 % Enable speed dependent losses
        Copper_Losses_Vary_With_Temperature (1,1) int32 % Enable copper loss variation with temperature
        RotorCopperLossesVaryWithTemp (1,1) int32 % Enable rotor cage losses variation with temperature
        StatorIronStrayLoadLossesVaryWithTemp (1,1) int32 % Enable stator iron stray load losses variation with temperature
        RotorIronStrayLoadLossesVaryWithTemp (1,1) int32 % Enable rotor iron stray load losses variation with temperature
        StatorCopperStrayLoadLossesVaryWithTemp (1,1) int32 % Enable stator copper stray load losses variation with temperature
        RotorCopperStrayLoadLossesVaryWithTemp (1,1) int32 % Enable rotor copper stray load losses variation with temperature
        ...
        Ventilated_Enable (1,1) int32 % Enable "Ventilated" cooling system
        Ventilated_FlowRate_m3ps (1,1) double % Flow rate (Ventilated) [m3/s]
        Ventilated_InletTemperature_degC (1,1) double % Inlet temperature (Ventilated) [degC]
        HousingWaterJacket_Enable (1,1) int32 % Enable "Housing Water Jacket" cooling system
        HousingWaterJacket_FlowRate_m3ps (1,1) double % Flow rate (Housing Water Jacket) [m3/s]
        HousingWaterJacket_InletTemperature_degC (1,1) double % Inlet temperature (Housing Water Jacket) [degC]
        ShaftSpiralGroove_Enable (1,1) int32  % Enable "Shaft Spiral Groove" cooling system
        ShaftSpiralGroove_FlowRate_m3ps (1,1) double % Flow rate (Shaft Spiral Groove) [m3/s]
        ShaftSpiralGroove_InletTemperature_degC (1,1) double % Inlet temperature (Shaft Spiral Groove) [degC]
        WetRotor_Enable (1,1) int32 % Enable "Wet Rotor" cooling system
        WetRotor_FlowRate_m3ps (1,1) double % Flow rate (Wet Rotor) [m3/s]
        WetRotor_InletTemperature_degC (1,1) double % Inlet temperature (Wet Rotor) [degC]
        SprayCooling_Enable (1,1) int32 % Enable "Spray Cooling" cooling system
        SprayCooling_FlowRate_m3ps (1,1) double % Flow rate (Spray Cooling) [m3/s]
        SprayCooling_InletTemperature_degC (1,1) double % Inlet temperature (Spray Cooling) [degC]
        RotorWaterJacket_Enable (1,1) int32 % Enable "Rotor Water Jacket" cooling system
        RotorWaterJacket_FlowRate_m3ps (1,1) double % Flow rate (Rotor Water Jacket) [m3/s]
        RotorWaterJacket_InletTemperature_degC (1,1) double % Inlet temperature (Rotor Water Jacket) [degC]
        SlotWaterJacket_Enable (1,1) int32 % Enable "Slot Water Jacket" cooling system
        SlotWaterJacket_FlowRate_m3ps (1,1) double % Flow rate (Slot Water Jacket) [m3/s]
        SlotWaterJacket_InletTemperature_degC (1,1) double % Inlet temperature (Slot Water Jacket) [degC]
        ...
        ThermalSteadyOrTransientChoice (1,1) int32 % Choice for steady-state or transient calculation
        TransientOption (1,1) int32 % Transient option
        StopTime (1,1) double % Transient stop time [s]
        StopTimeChoice (1,1) int32 % Choice for transient stop condition
        TransientLossesOrTorqueChoice (1,1) int32 % Choice for transient losses (specify loss or specify torque)
        NumTimeSteps (1,1) int32 % Number of time steps for transient calculation
        TransientTorqueValue (1,1) double % Load torque for transient calculation
        InitialTemperatureOption (1,1) int32 % Initial temperature option for transient calculation
    end
    properties(GetAccess=protected)
        Ventilated_Enable2 (1,1) int32
        Shaft_Speed_RPM_Thermal (1,1) double
    end
    % Losses
    properties(GetAccess=protected)                
        Armature_Copper_Loss (1,1) double % Stator copper loss [W]
        Rotor_Copper_Loss (1,1) double % Rotor copper loss [W]
        Stator_Iron_Loss_Back_Iron (1,1) double % Stator back iron loss [W]
        Stator_Iron_Loss_Tooth (1,1) double % Stator tooth iron loss [W]
        Stator_Iron_Stray_Load_Loss (1,1) double % Stator iron stray loss [W]
        Rotor_Iron_Stray_Load_Loss (1,1) double % Rotor iron stray loss [W]
        Stator_Copper_Stray_Load_Loss (1,1) double % Stator copper stray loss [W]
        Rotor_Copper_Stray_Load_Loss (1,1) double % Rotor copper stray loss [W]
        Magnet_Iron_Loss (1,1) double % Magnet loss [W]
        Rotor_Iron_Loss_Embedded_Magnet_Pole (1,1) double % Rotor pole iron loss [W]
        Rotor_Iron_Loss_Back_Iron (1,1) double % Rotor back iron loss [W]
        Rotor_Iron_Loss_Tooth (1,1) double % Rotor tooth iron loss [W]
        Friction_Loss_F (1,1) double % Front bearing friction loss [W]
        Friction_Loss_R (1,1) double % Rear bearing friction loss [W]
        Windage_Loss (1,1) double % Windage loss [W]
        Windage_Loss_Ext_Fan (1,1) double % External fan windage loss [W]
    end
    properties(Access=public)
        Stray_Loss_Stator_Iron_Proportion (1,1) double % Proportion of iron stray loss in stator
    end
    properties(Constant, Access=public)
        % LossTypes
        LossTypes = {'Armature_Copper_Loss'; ...
                     'Rotor_Copper_Loss'; ...
                     'Stator_Iron_Loss_Back_Iron'; ...
                     'Stator_Iron_Loss_Tooth'; ...
                     'Stator_Iron_Stray_Load_Loss'; ...
                     'Rotor_Iron_Stray_Load_Loss'; ...
                     'Stator_Copper_Stray_Load_Loss'; ...
                     'Rotor_Copper_Stray_Load_Loss'; ...
                     'Magnet_Iron_Loss'; ...
                     'Rotor_Iron_Loss_Embedded_Magnet_Pole'; ...
                     'Rotor_Iron_Loss_Back_Iron'; ...
                     'Rotor_Iron_Loss_Tooth'; ...
                     'Friction_Loss_F'; ...
                     'Friction_Loss_R'; ...
                     'Windage_Loss'; ...
                     'Windage_Loss_Ext_Fan'} 
    end
    properties(Access=public)
        LossValues (16,1) double % Loss values corresponding to LossTypes
    end

    % Overload the constructor
    methods
        function obj = BasicInterface(motFullFile)
            % BasicInterface constructor

            obj = obj@mcadinterface.AbstractInterface(motFullFile);

            % Special parameters
            for idxLoss = 1:length(obj.LossTypes)
                thisLossName = obj.LossTypes{idxLoss};
                thisLossValue = obj.getParameter(thisLossName);
                obj.LossValues(idxLoss) = thisLossValue;
            end

            % SelfVentilated OR Through_Ventilated are the same setting. Use Through_Ventilated only
            enableSelfVent = obj.getParameter('Ventilated_Enable2');
            if enableSelfVent
                obj.Ventilated_Enable2 = 0;
                obj.Ventilated_Enable = 1;
            end

        end
    end

    % Basic set methods
    methods
        function set.Tambient_degC(obj, value)
            obj.setParameter('Tambient_degC', value);
            obj.Tambient_degC = value;
        end

        function set.EnableStatorTempCoeffRes(obj, value)
            obj.setParameter('EnableStatorTempCoeffRes', value);
            obj.EnableStatorTempCoeffRes = value;
        end

        function set.StatorTempCoeffRes(obj, value)
            obj.setParameter('StatorTempCoeffRes', value);
            obj.StatorTempCoeffRes = value;
        end

        function set.TrefStator_degC(obj, value)
            obj.setParameter('TrefStator_degC', value);
            obj.TrefStator_degC = value;
        end

        function set.EnableRotorTempCoeffRes(obj, value)
            obj.setParameter('EnableRotorTempCoeffRes', value);
            obj.EnableRotorTempCoeffRes = value;
        end

        function set.RotorTempCoeffRes(obj, value)
            obj.setParameter('RotorTempCoeffRes', value);
            obj.RotorTempCoeffRes = value;
        end

        function set.TrefRotor_degC(obj, value)
            obj.setParameter('TrefRotor_degC', value);
            obj.TrefRotor_degC = value;
        end

        function set.MagneticThermalCouplingChoice(obj, value)
            obj.setParameter('MagneticThermalCouplingChoice', value);
            obj.MagneticThermalCouplingChoice = value;
        end

        function set.LabThermalCouplingChoice(obj, value)
            obj.setParameter('LabThermalCouplingChoice', value);
            obj.LabThermalCouplingChoice = value;
        end

        function set.Armature_Copper_Loss(obj, value)
            obj.setParameter('Armature_Copper_Loss', value);
            obj.Armature_Copper_Loss = value;
        end

        function set.Rotor_Copper_Loss(obj, value)
            obj.setParameter('Rotor_Copper_Loss', value);
            obj.Rotor_Copper_Loss = value;
        end

        function set.Stator_Iron_Loss_Back_Iron(obj, value)
            obj.setParameter('Stator_Iron_Loss_Back_Iron', value);
            obj.Stator_Iron_Loss_Back_Iron = value;
        end

        function set.Stator_Iron_Loss_Tooth(obj, value)
            obj.setParameter('Stator_Iron_Loss_Tooth', value);
            obj.Stator_Iron_Loss_Tooth = value;
        end

        function set.Stator_Iron_Stray_Load_Loss(obj, value)
            obj.setParameter('Stator_Iron_Stray_Load_Loss', value);
            obj.Stator_Iron_Stray_Load_Loss = value;
        end

        function set.Rotor_Iron_Stray_Load_Loss(obj, value)
            obj.setParameter('Rotor_Iron_Stray_Load_Loss', value);
            obj.Rotor_Iron_Stray_Load_Loss = value;
        end

        function set.Stator_Copper_Stray_Load_Loss(obj, value)
            obj.setParameter('Stator_Copper_Stray_Load_Loss', value);
            obj.Stator_Copper_Stray_Load_Loss = value;
        end

        function set.Rotor_Copper_Stray_Load_Loss(obj, value)
            obj.setParameter('Rotor_Copper_Stray_Load_Loss', value);
            obj.Rotor_Copper_Stray_Load_Loss = value;
        end

        function set.Magnet_Iron_Loss(obj, value)
            obj.setParameter('Magnet_Iron_Loss', value);
            obj.Magnet_Iron_Loss = value;
        end

        function set.Rotor_Iron_Loss_Embedded_Magnet_Pole(obj, value)
            obj.setParameter('Rotor_Iron_Loss_Embedded_Magnet_Pole', value);
            obj.Rotor_Iron_Loss_Embedded_Magnet_Pole = value;
        end

        function set.Rotor_Iron_Loss_Back_Iron(obj, value)
            obj.setParameter('Rotor_Iron_Loss_Back_Iron', value);
            obj.Rotor_Iron_Loss_Back_Iron = value;
        end

        function set.Friction_Loss_F(obj, value)
            obj.setParameter('Friction_Loss_F', value);
            obj.Friction_Loss_F = value;
        end

        function set.Friction_Loss_R(obj, value)
            obj.setParameter('Friction_Loss_R', value);
            obj.Friction_Loss_R = value;
        end

        function set.Windage_Loss(obj, value)
            obj.setParameter('Windage_Loss', value);
            obj.Windage_Loss = value;
        end

        function set.Windage_Loss_Ext_Fan(obj, value)
            obj.setParameter('Windage_Loss_Ext_Fan', value);
            obj.Windage_Loss_Ext_Fan = value;
        end

        function set.Stray_Loss_Stator_Iron_Proportion(obj, value)
            obj.setParameter('Stray_Loss_Stator_Iron_Proportion', value);
            obj.Stray_Loss_Stator_Iron_Proportion = value;
        end       

        function set.Shaft_Torque_Nm(obj, value)
            obj.setParameter('Shaft_Torque_Nm', value);
            obj.Shaft_Torque_Nm = value;
        end

        % --------------
        function set.Loss_Function_Speed(obj, value)
            obj.setParameter('Loss_Function_Speed', value);
            obj.Loss_Function_Speed = value;
        end

        function set.Copper_Losses_Vary_With_Temperature(obj, value)
            obj.setParameter('Copper_Losses_Vary_With_Temperature', value);
            obj.Copper_Losses_Vary_With_Temperature = value;
        end

        function set.RotorCopperLossesVaryWithTemp(obj, value)
            obj.setParameter('RotorCopperLossesVaryWithTemp', value);
            obj.RotorCopperLossesVaryWithTemp = value;
        end

        function set.StatorIronStrayLoadLossesVaryWithTemp(obj, value)
            obj.setParameter('StatorIronStrayLoadLossesVaryWithTemp', value);
            obj.StatorIronStrayLoadLossesVaryWithTemp = value;
        end

        function set.RotorIronStrayLoadLossesVaryWithTemp(obj, value)
            obj.setParameter('RotorIronStrayLoadLossesVaryWithTemp', value);
            obj.RotorIronStrayLoadLossesVaryWithTemp = value;
        end

        function set.StatorCopperStrayLoadLossesVaryWithTemp(obj, value)
            obj.setParameter('StatorCopperStrayLoadLossesVaryWithTemp', value);
            obj.StatorCopperStrayLoadLossesVaryWithTemp = value;
        end

        function set.RotorCopperStrayLoadLossesVaryWithTemp(obj, value)
            obj.setParameter('RotorCopperStrayLoadLossesVaryWithTemp', value);
            obj.RotorCopperStrayLoadLossesVaryWithTemp = value;
        end


        % ---------------

        function set.Ventilated_Enable(obj, value)
            obj.setParameter('Ventilated_Enable', value); 
            obj.Ventilated_Enable = value;
        end

        function set.Ventilated_FlowRate_m3ps(obj, value)
            obj.setParameter('Ventilated_FlowRate_m3ps', value);
            obj.Ventilated_FlowRate_m3ps = value;
        end

        function set.Ventilated_InletTemperature_degC(obj, value)
            obj.setParameter('Ventilated_InletTemperature_degC', value);
            obj.Ventilated_InletTemperature_degC = value;
        end

        function set.HousingWaterJacket_Enable(obj, value)
            obj.setParameter('HousingWaterJacket_Enable', value);
            obj.HousingWaterJacket_Enable = value;
        end

        function set.HousingWaterJacket_FlowRate_m3ps(obj, value)
            obj.setParameter('HousingWaterJacket_FlowRate_m3ps', value);
            obj.HousingWaterJacket_FlowRate_m3ps = value;
        end

        function set.HousingWaterJacket_InletTemperature_degC(obj, value)
            obj.setParameter('HousingWaterJacket_InletTemperature_degC', value);
            obj.HousingWaterJacket_InletTemperature_degC = value;
        end

        function set.ShaftSpiralGroove_Enable(obj, value)
            obj.setParameter('ShaftSpiralGroove_Enable', value);
            obj.ShaftSpiralGroove_Enable = value;
        end

        function set.ShaftSpiralGroove_FlowRate_m3ps(obj, value)
            obj.setParameter('ShaftSpiralGroove_FlowRate_m3ps', value);
            obj.ShaftSpiralGroove_FlowRate_m3ps = value;
        end

        function set.ShaftSpiralGroove_InletTemperature_degC(obj, value)
            obj.setParameter('ShaftSpiralGroove_InletTemperature_degC', value);
            obj.ShaftSpiralGroove_InletTemperature_degC = value;
        end

        function set.WetRotor_Enable(obj, value)
            obj.setParameter('WetRotor_Enable', value);
            obj.WetRotor_Enable = value;
        end

        function set.WetRotor_FlowRate_m3ps(obj, value)
            obj.setParameter('WetRotor_FlowRate_m3ps', value);
            obj.WetRotor_FlowRate_m3ps = value;
        end

        function set.WetRotor_InletTemperature_degC(obj, value)
            obj.setParameter('WetRotor_InletTemperature_degC', value);
            obj.WetRotor_InletTemperature_degC = value;
        end

        function set.SprayCooling_Enable(obj, value)
            obj.setParameter('SprayCooling_Enable', value);
            obj.SprayCooling_Enable = value;
        end

        function set.SprayCooling_FlowRate_m3ps(obj, value)
            obj.setParameter('SprayCooling_FlowRate_m3ps', value);
            obj.SprayCooling_FlowRate_m3ps = value;
        end

        function set.SprayCooling_InletTemperature_degC(obj, value)
            obj.setParameter('SprayCooling_InletTemperature_degC', value);
            obj.SprayCooling_InletTemperature_degC = value;
        end

        function set.RotorWaterJacket_Enable(obj, value)
            obj.setParameter('RotorWaterJacket_Enable', value);
            obj.RotorWaterJacket_Enable = value;
        end

        function set.RotorWaterJacket_FlowRate_m3ps(obj, value)
            obj.setParameter('RotorWaterJacket_FlowRate_m3ps', value);
            obj.RotorWaterJacket_FlowRate_m3ps = value;
        end

        function set.RotorWaterJacket_InletTemperature_degC(obj, value)
            obj.setParameter('RotorWaterJacket_InletTemperature_degC', value);
            obj.RotorWaterJacket_InletTemperature_degC = value;
        end

        function set.SlotWaterJacket_Enable(obj, value)
            obj.setParameter('SlotWaterJacket_Enable', value);
            obj.SlotWaterJacket_Enable = value;
        end

        function set.SlotWaterJacket_FlowRate_m3ps(obj, value)
            obj.setParameter('SlotWaterJacket_FlowRate_m3ps', value);
            obj.SlotWaterJacket_FlowRate_m3ps = value;
        end

        function set.SlotWaterJacket_InletTemperature_degC(obj, value)
            obj.setParameter('SlotWaterJacket_InletTemperature_degC', value);
            obj.SlotWaterJacket_InletTemperature_degC = value;
        end

        function set.ThermalSteadyOrTransientChoice(obj, value)
            obj.setParameter('ThermalSteadyOrTransientChoice', value);
            obj.ThermalSteadyOrTransientChoice = value;
        end

        function set.TransientOption(obj, value)
            obj.setParameter('TransientOption', value);
            obj.TransientOption = value;
        end

        function set.StopTime(obj, value)
            obj.setParameter('StopTime', value);
            obj.StopTime = value;
        end

        function set.StopTimeChoice(obj, value)
            obj.setParameter('StopTimeChoice', value);
            obj.StopTimeChoice = value;
        end

        function set.TransientLossesOrTorqueChoice(obj, value)
            obj.setParameter('TransientLossesOrTorqueChoice', value);
            obj.TransientLossesOrTorqueChoice = value;
        end

        function set.NumTimeSteps(obj, value)
            obj.setParameter('NumTimeSteps', value);
            obj.NumTimeSteps = value;
        end

        function set.TransientTorqueValue(obj, value)
            obj.setParameter('TransientTorqueValue', value);
            obj.TransientTorqueValue = value;
        end

        function set.InitialTemperatureOption(obj, value)
            obj.setParameter('InitialTemperatureOption', value);
            obj.InitialTemperatureOption = value;
        end        
    
    end

    % Composite set methods
    methods

        function set.Shaft_Speed_RPM(obj, value)
            obj.setParameter('Shaft_Speed_RPM', value);
            obj.Shaft_Speed_RPM = value;
            obj.setParameter('Shaft_Speed_RPM_Thermal', value);
        end

        function set.LossValues(obj, valuesVec)
            arguments
                obj mcadinterface.BasicInterface
                valuesVec (:,1) {mustBeNumeric,mustBeGreaterThanOrEqual(valuesVec, 0)}
            end
            for idxLoss = 1:length(obj.LossTypes)
                thisLossName = obj.LossTypes{idxLoss};
                obj.setParameter(thisLossName, valuesVec(idxLoss))
            end
            obj.LossValues = valuesVec;
        end

    end
end


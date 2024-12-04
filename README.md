# **Import a Motor-CAD Thermal Model into Simulink and Simscape**

[![View Import a Motor-CAD Thermal Model into Simulink and Simscape on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/120598-import-a-motor-cad-thermal-model-into-simulink-and-simscape) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=mathworks/import-motorcad-thermal-simulink)

## Overview
You can use this project to import an Ansys™ Motor-CAD™ motor model into Simulink®. 
This version of the project is intended for versions of Motor-CAD that use the Visual Basic (ActiveX) interface.
For newer Motor-CAD versions, such as 2023R1 or newer, which use the PyMotorCAD interface, use a newer version. See all versions in [GitHub Releases](https://github.com/mathworks/import-motorcad-thermal-simulink/releases).

You can use the automatically-generated Simulink model to:
 * Predict the transient temperature of the motor elements under dynamic operating points and diverse cooling scenarios. 
 * Run simulations faster than real time.
 * Integrate the motor in a system-level model using Simscape™.
 * Analize the performance of your system design in a holistic way.

## 1. Import a Motor-CAD model into MATLAB ##
You can use the Object-Oriented MATLAB® - Motor-CAD interface (included in this repository) to import a Motor-CAD model into MATLAB and easily modify motor properties and run Motor-CAD calculations.

For more information, run >> doc mcadinterface.ThermalInterface

## 2. Generate a Simulink Reduced-Order Thermal Model (SROTM) ##
Run the *GenerateSimulinkThermalModel.mlx* live script to automatically generate the SROTM for an induction motor (IM) and a permanent magnet synchronous motor (PMSM).

The IM has one active cooling system (housing water jacket). It is based on the *e5_eMobility_IM* Motor-CAD template.

The PMSM has two active cooling systems (housing water jacket and through-ventilation). It is based on the *e8_eMobility_IPM* Motor-CAD template.

![](/images/e8_IPMSM_HWJandVent_ROM_screenshot.PNG)

## 3. Validate the SROTM ##
Run the *ValidateSimulinkThermalModel.mlx* live script to compare the SROTM simulation results with the Motor-CAD simulation results at different operating points.

![](/images/validateSROTM_screenshot.PNG)

## 4. Use the SROTM in Simscape Models ##
Run the *UseSimulinkThermalModel.mlx* live script to simulate a Simscape system-level vehicle model that integrates the Simulink reduced-order thermal model of the motor.

![](/images/SimscapeSystemLevel_screenshot.PNG)

## Setup 
Open the project file *ImportMotorCADThermalModel.prj* to get started.
- Run *GenerateSimulinkThermalModel.mlx* (requires Motor-CAD)
- Run *ValidateSimulinkThermalModel.mlx* (requires Motor-CAD)
- Run *UseSimulinkThermalModel.mlx* 

### MathWorks Products (http://www.mathworks.com)
Requires MATLAB® release R2021b or newer.
- [Simulink™](https://www.mathworks.com/products/simulink.html) - required.
- [Simscape™](https://www.mathworks.com/products/simscape.html) - optional.
- [Simscape™ Electrical™](https://www.mathworks.com/products/simscape-electrical.html) - optional.
- [Simscape™ Fluids™](https://www.mathworks.com/products/simscape-fluids.html) - optional.
- [Control System Toolbox™](https://www.mathworks.com/products/control.html) - required (LPV System block).
- [Motor Control Blockset™](https://www.mathworks.com/products/motor-control.html) - optional (Drive Cycle Source block).

### Getting Started 
To learn more about modeling and simulation with Simscape™, please visit:
* [Simscape Getting Started Resources](https://www.mathworks.com/solutions/physical-modeling/resources.html)

## License
The license is available in the License file within this repository.

## Community Support
[MATLAB Central](https://www.mathworks.com/matlabcentral)

Copyright 2022 The MathWorks, Inc.

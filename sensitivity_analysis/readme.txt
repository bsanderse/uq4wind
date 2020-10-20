%% Sensitivity analysis and calibration for ECN AeroModule

In order to run the code, you need 
- UQLab and a license for UQLab. Set your UQLab path in config.m. Note that version 1.3 is provided standard with this repository (but not the license).
- ECN AeroModule and a license for it (contact Koen Boorsma at TNO). ECNAero.exe should exist on your path (set this in environment variables on WindowS)

For the AeroModule runs, a folder that defines the test case should be present in the cases/ folder
in this folder, let's say XXX, one needs the following files:
- initialize.m or initialize_calibration.m (depending whether Sensitivity or Calibration is performed)
- XXX.m containing the definition of uncertain parameters
- postprocessing.m or postprocessing_calibration.m
- XXX_readoutput.m which reads the AeroModule output and selects the QoI
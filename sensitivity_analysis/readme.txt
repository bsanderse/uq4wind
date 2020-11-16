%% Sensitivity analysis and calibration for ECN AeroModule

In order to run the code, you need 
- UQLab and a license for UQLab. Set your UQLab path in config.m. Note that version 1.3 is provided standard with this repository (but not the license).
- ECN AeroModule and a license for it (contact Koen Boorsma at TNO). ECNAero.exe should exist on your path (set this in environment variables on Windows). 
You can use the function findAeroModulePath() (in the folder Other/) to check where ECNAero.exe is located on your path.

For the AeroModule runs, a folder that defines the test case should be present in the cases/ folder.
In this folder, let's say XXX, one needs the following files:
- an initialization file: initialize.m or initialize_calibration.m (depending whether Sensitivity or Calibration is performed)
- XXX.m containing the definition of uncertain parameters
- XXX_readoutput.m which reads the AeroModule output and selects the quantity of interest (QoI)
- a postprocessing file, postprocessing.m or postprocessing_calibration.m

You can then run the code as follows:
- For sensitivity analysis, run testSensitivity.m, by specifying the casename in testSensitivity.m
- For simple calibration, use testCalibration.m, by specifying the casename in testCalibration.m
- For advanced calibration, with multiple operating conditions, use Calibration_multipleruns.m, and specify the casename in Calibration_multipleruns.m

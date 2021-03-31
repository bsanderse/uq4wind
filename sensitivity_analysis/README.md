# Sensitivity analysis and calibration using UQ4WIND and UQLab

In order to use UQ4WIND you need:
- UQLab and a license for UQLab. Set your UQLab path in config.m. Note that version 1.3 is provided standard with this repository (but not the license). The latest UQLab version is 1.4 and should work as well, although UQ4WIND has been tested mostly with v1.3.

Setting up a test case is done by creating a folder inside the cases/ folder.
In this folder one needs the following files:
- An initialization file: initialize.m or initialize_calibration.m (depending whether Sensitivity or Calibration is performed)
- A postprocessing file, postprocessing.m or postprocessing_calibration.m

You can then run the code as follows:
- For sensitivity analysis, go into SensitivityAnalysis.m and specify the casename. Then run SensitivityAnalysis.m.
- For simple calibration, go into Calibration.m and specify the casename. Then run Calibration.m.
- For advanced calibration, with multiple operating conditions, use Calibration_multipleruns.m

## Examples
Simple examples for sensitivity analysis that should run 'out of the box' are
- airfoil_lift
- Cinf
- Ishigami
- linear_portfolio

## AeroModule runs
For the AeroModule runs, you additionally need:
- ECN AeroModule and a license for it (contact Koen Boorsma at TNO). ECNAero.exe should exist on your path (set this in environment variables on Windows). 
You can use the function findAeroModulePath() (in the folder Other/) to check where ECNAero.exe is located on your path.

In addition to the initialization and postprocessing file, for the AeroModule runs one needs (assuming the casename is XXX):
- XXX.m containing the definition of uncertain parameters
- XXX_readoutput.m which reads the AeroModule output and selects the quantity of interest (QoI)

## Workflow diagram
<img src="workflow_windtrue.png" width="800" height="450">

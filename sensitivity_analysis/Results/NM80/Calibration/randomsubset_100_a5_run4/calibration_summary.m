%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            100

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            100

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            100

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            100

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       4

%------------------- Solver
   Solution method:                               MCMC

   Algorithm:                                     AIES
   Duration (HH:MM:SS):                           00:14:31
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.22   | 0.016   | (-0.24 - -0.19)     | Model       |
| CL        | -0.15   | 0.019   | (-0.19 - -0.12)     | Model       |
| CL        | -0.15   | 0.056   | (-0.23 - -0.045)    | Model       |
| CL        | -0.2    | 0.0088  | (-0.21 - -0.19)     | Model       |
| Sigma2-1  | 3.3e+03 | 4.7e+02 | (2.6e+03 - 4.1e+03) | Discrepancy |
| Sigma2-2  | 6.5e+03 | 9.5e+02 | (5.1e+03 - 8.2e+03) | Discrepancy |
| Sigma2-3  | 1.2e+04 | 1.9e+03 | (9.7e+03 - 1.6e+04) | Discrepancy |
| Sigma2-4  | 1e+04   | 1.4e+03 | (8e+03 - 1.3e+04)   | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.22   | Model          |
| CL        | -0.15   | Model          |
| CL        | -0.15   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 3.3e+03 | Discrepancy    |
| Sigma2-2  | 6.5e+03 | Discrepancy    |
| Sigma2-3  | 1.2e+04 | Discrepancy    |
| Sigma2-4  | 1e+04   | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
----------------------------------------------
|    |  CL        CL       CL        CL      |
----------------------------------------------
| CL |  1         -0.14    -0.016    -0.065  |
| CL |  -0.14     1        -0.24     0.15    |
| CL |  -0.016    -0.24    1         -0.47   |
| CL |  -0.065    0.15     -0.47     1       |
----------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.0022     0.078       -0.0095   |
| Sigma2-2 |  -0.0022     1           -0.023      0.029     |
| Sigma2-3 |  0.078       -0.023      1           -0.036    |
| Sigma2-4 |  -0.0095     0.029       -0.036      1         |
-------------------------------------------------------------



The MCMC simulation has converged
Warning: unperturbed AeroModule results are obtained from the surrogate model 
> In postProcessing_calibration (line 68)
  In run (line 91)
  In testCalibration (line 137) 
note: no specialist input file
note: no specialist input file
Warning: Could not get change notification handle for
\\mac\home\Dropbox\work\Programming\UQ\windtrue\sensitivity_analysis\Results\NM80\Calibration.
Type 'help changeNotification' for more info. 
for i=1:12
figure(i);savefig(['fig' num2str(i)]);
end

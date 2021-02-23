The estimation of PCE coefficients stopped at polynomial degree 4 and qNorm 0.75 for output variable 1
Final LOO error estimate: 3.068386e-06
---                 Calculation finished!                               ---
---   Calculating the PCE coefficients with least-squares.   ---
The estimation of PCE coefficients stopped at polynomial degree 4 and qNorm 0.75 for output variable 2
Final LOO error estimate: 4.750837e-05
---                 Calculation finished!                               ---
---   Calculating the PCE coefficients with least-squares.   ---
The estimation of PCE coefficients stopped at polynomial degree 4 and qNorm 0.75 for output variable 3
Final LOO error estimate: 6.814234e-05
---                 Calculation finished!                               ---
---   Calculating the PCE coefficients with least-squares.   ---
The estimation of PCE coefficients stopped at polynomial degree 3 and qNorm 0.75 for output variable 4
Final LOO error estimate: 2.183688e-05
---                 Calculation finished!                               ---
performing Bayesian analysis

Starting AIES...

|##############################| 100.00%

Finished AIES!

%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            1000

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            1000

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            1000

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            1000

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
   Duration (HH:MM:SS):                           00:49:53
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.23   | 0.0056  | (-0.24 - -0.23)     | Model       |
| CL        | -0.17   | 0.0062  | (-0.18 - -0.16)     | Model       |
| CL        | -0.21   | 0.0037  | (-0.22 - -0.21)     | Model       |
| CL        | -0.2    | 0.0034  | (-0.21 - -0.2)      | Model       |
| Sigma2-1  | 4.1e+03 | 1.8e+02 | (3.8e+03 - 4.4e+03) | Discrepancy |
| Sigma2-2  | 8e+03   | 3.6e+02 | (7.4e+03 - 8.6e+03) | Discrepancy |
| Sigma2-3  | 1.4e+04 | 6.5e+02 | (1.3e+04 - 1.5e+04) | Discrepancy |
| Sigma2-4  | 1.2e+04 | 5.1e+02 | (1.1e+04 - 1.2e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.23   | Model          |
| CL        | -0.17   | Model          |
| CL        | -0.21   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 4.1e+03 | Discrepancy    |
| Sigma2-2  | 8e+03   | Discrepancy    |
| Sigma2-3  | 1.4e+04 | Discrepancy    |
| Sigma2-4  | 1.2e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
----------------------------------------------
|    |  CL       CL        CL        CL      |
----------------------------------------------
| CL |  1        -0.17     0.036     0.014   |
| CL |  -0.17    1         -0.038    -0.036  |
| CL |  0.036    -0.038    1         -0.17   |
| CL |  0.014    -0.036    -0.17     1       |
----------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           0.014       -0.015      -0.0024   |
| Sigma2-2 |  0.014       1           -0.021      0.01      |
| Sigma2-3 |  -0.015      -0.021      1           0.0023    |
| Sigma2-4 |  -0.0024     0.01        0.0023      1         |
-------------------------------------------------------------



The MCMC simulation has converged
Warning: unperturbed AeroModule results are obtained from the surrogate model 
> In postProcessing_calibration (line 68)
  In run (line 91)
  In testCalibration (line 137) 
note: no specialist input file
note: no specialist input file
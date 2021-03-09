%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            10

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            10

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            10

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            10

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
   Duration (HH:MM:SS):                           00:12:33
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.21   | 0.055   | (-0.27 - -0.098)    | Model       |
| CL        | -0.15   | 0.037   | (-0.2 - -0.086)     | Model       |
| CL        | -0.22   | 0.02    | (-0.25 - -0.19)     | Model       |
| CL        | -0.2    | 0.022   | (-0.24 - -0.16)     | Model       |
| Sigma2-1  | 2.2e+03 | 2.9e+03 | (4.8e+02 - 7.4e+03) | Discrepancy |
| Sigma2-2  | 2.2e+03 | 2.7e+03 | (5e+02 - 6.9e+03)   | Discrepancy |
| Sigma2-3  | 3.6e+03 | 2.8e+03 | (1.3e+03 - 8.4e+03) | Discrepancy |
| Sigma2-4  | 5.6e+03 | 5.9e+03 | (1.4e+03 - 1.8e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.21   | Model          |
| CL        | -0.15   | Model          |
| CL        | -0.22   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 2.2e+03 | Discrepancy    |
| Sigma2-2  | 2.2e+03 | Discrepancy    |
| Sigma2-3  | 3.6e+03 | Discrepancy    |
| Sigma2-4  | 5.6e+03 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-------------------------------------------------
|    |  CL        CL         CL        CL       |
-------------------------------------------------
| CL |  1         -0.11      0.07      -0.058   |
| CL |  -0.11     1          -0.021    -0.0077  |
| CL |  0.07      -0.021     1         -0.14    |
| CL |  -0.058    -0.0077    -0.14     1        |
-------------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.11       -0.015      -0.04     |
| Sigma2-2 |  -0.11       1           -0.059      -0.067    |
| Sigma2-3 |  -0.015      -0.059      1           -0.016    |
| Sigma2-4 |  -0.04       -0.067      -0.016      1         |
-------------------------------------------------------------




R_hat_full =

    2.3996
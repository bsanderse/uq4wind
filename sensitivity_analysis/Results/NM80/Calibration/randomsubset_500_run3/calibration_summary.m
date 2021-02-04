%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            500

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            500

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            500

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            500

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
   Duration (HH:MM:SS):                           00:42:34
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.22   | 0.0086  | (-0.24 - -0.21)     | Model       |
| CL        | -0.17   | 0.0094  | (-0.18 - -0.15)     | Model       |
| CL        | -0.13   | 0.03    | (-0.18 - -0.079)    | Model       |
| CL        | -0.2    | 0.0038  | (-0.2 - -0.19)      | Model       |
| Sigma2-1  | 3.9e+03 | 2.5e+02 | (3.5e+03 - 4.3e+03) | Discrepancy |
| Sigma2-2  | 7.4e+03 | 4.4e+02 | (6.7e+03 - 8.2e+03) | Discrepancy |
| Sigma2-3  | 1.4e+04 | 9.1e+02 | (1.3e+04 - 1.6e+04) | Discrepancy |
| Sigma2-4  | 1.1e+04 | 7e+02   | (9.8e+03 - 1.2e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.22   | Model          |
| CL        | -0.17   | Model          |
| CL        | -0.13   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 3.9e+03 | Discrepancy    |
| Sigma2-2  | 7.4e+03 | Discrepancy    |
| Sigma2-3  | 1.4e+04 | Discrepancy    |
| Sigma2-4  | 1.1e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
---------------------------------------------
|    |  CL        CL       CL       CL      |
---------------------------------------------
| CL |  1         -0.12    0.029    -0.048  |
| CL |  -0.12     1        -0.17    0.12    |
| CL |  0.029     -0.17    1        -0.47   |
| CL |  -0.048    0.12     -0.47    1       |
---------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           0.007       0.0069      0.035     |
| Sigma2-2 |  0.007       1           0.021       -0.02     |
| Sigma2-3 |  0.0069      0.021       1           0.036     |
| Sigma2-4 |  0.035       -0.02       0.036       1         |
-------------------------------------------------------------



The MCMC simulation has converged
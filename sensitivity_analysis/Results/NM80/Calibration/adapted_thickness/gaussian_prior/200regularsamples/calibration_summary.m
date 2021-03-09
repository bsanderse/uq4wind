%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            204

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            204

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            204

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            204

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
   Duration (HH:MM:SS):                           00:27:05
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.24   | 0.013   | (-0.26 - -0.21)     | Model       |
| CL        | -0.18   | 0.014   | (-0.2 - -0.15)      | Model       |
| CL        | -0.23   | 0.0091  | (-0.24 - -0.21)     | Model       |
| CL        | -0.21   | 0.0075  | (-0.22 - -0.19)     | Model       |
| Sigma2-1  | 4e+03   | 4e+02   | (3.4e+03 - 4.8e+03) | Discrepancy |
| Sigma2-2  | 8.2e+03 | 8.3e+02 | (7e+03 - 9.7e+03)   | Discrepancy |
| Sigma2-3  | 1.5e+04 | 1.5e+03 | (1.3e+04 - 1.8e+04) | Discrepancy |
| Sigma2-4  | 1.3e+04 | 1.3e+03 | (1.1e+04 - 1.5e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.24   | Model          |
| CL        | -0.18   | Model          |
| CL        | -0.23   | Model          |
| CL        | -0.21   | Model          |
| Sigma2-1  | 4e+03   | Discrepancy    |
| Sigma2-2  | 8.2e+03 | Discrepancy    |
| Sigma2-3  | 1.5e+04 | Discrepancy    |
| Sigma2-4  | 1.3e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-------------------------------------------------
|    |  CL         CL        CL         CL      |
-------------------------------------------------
| CL |  1          -0.086    -0.0038    0.0065  |
| CL |  -0.086     1         -0.023     -0.012  |
| CL |  -0.0038    -0.023    1          -0.12   |
| CL |  0.0065     -0.012    -0.12      1       |
-------------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           0.082       0.025       0.059     |
| Sigma2-2 |  0.082       1           -0.013      0.033     |
| Sigma2-3 |  0.025       -0.013      1           0.028     |
| Sigma2-4 |  0.059       0.033       0.028       1         |
-------------------------------------------------------------



The MCMC simulation has converged
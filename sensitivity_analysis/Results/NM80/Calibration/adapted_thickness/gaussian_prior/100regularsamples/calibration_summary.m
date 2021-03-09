%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            101

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            101

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            101

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            101

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
   Duration (HH:MM:SS):                           00:26:01
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.22   | 0.018   | (-0.25 - -0.19)     | Model       |
| CL        | -0.14   | 0.02    | (-0.17 - -0.1)      | Model       |
| CL        | -0.21   | 0.012   | (-0.23 - -0.19)     | Model       |
| CL        | -0.2    | 0.0099  | (-0.22 - -0.18)     | Model       |
| Sigma2-1  | 4.1e+03 | 6e+02   | (3.2e+03 - 5.2e+03) | Discrepancy |
| Sigma2-2  | 7.7e+03 | 1.2e+03 | (6e+03 - 9.8e+03)   | Discrepancy |
| Sigma2-3  | 1.2e+04 | 1.7e+03 | (9.7e+03 - 1.5e+04) | Discrepancy |
| Sigma2-4  | 1.1e+04 | 1.7e+03 | (8.9e+03 - 1.4e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.22   | Model          |
| CL        | -0.14   | Model          |
| CL        | -0.21   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 4.1e+03 | Discrepancy    |
| Sigma2-2  | 7.7e+03 | Discrepancy    |
| Sigma2-3  | 1.2e+04 | Discrepancy    |
| Sigma2-4  | 1.1e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-----------------------------------------------
|    |  CL        CL        CL        CL      |
-----------------------------------------------
| CL |  1         -0.15     -0.02     0.0038  |
| CL |  -0.15     1         -0.024    -0.046  |
| CL |  -0.02     -0.024    1         -0.16   |
| CL |  0.0038    -0.046    -0.16     1       |
-----------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.00081    -0.033      0.089     |
| Sigma2-2 |  -0.00081    1           -0.0012     -0.017    |
| Sigma2-3 |  -0.033      -0.0012     1           0.0061    |
| Sigma2-4 |  0.089       -0.017      0.0061      1         |
-------------------------------------------------------------



The MCMC simulation has converged
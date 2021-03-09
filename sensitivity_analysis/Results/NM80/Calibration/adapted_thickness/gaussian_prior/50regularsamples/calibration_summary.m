%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            51

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            51

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            51

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            51

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
   Duration (HH:MM:SS):                           00:16:04
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.2    | 0.027   | (-0.24 - -0.15)     | Model       |
| CL        | -0.16   | 0.026   | (-0.2 - -0.12)      | Model       |
| CL        | -0.21   | 0.018   | (-0.24 - -0.18)     | Model       |
| CL        | -0.2    | 0.014   | (-0.22 - -0.17)     | Model       |
| Sigma2-1  | 4e+03   | 8.5e+02 | (2.9e+03 - 5.6e+03) | Discrepancy |
| Sigma2-2  | 6.9e+03 | 1.4e+03 | (5e+03 - 9.6e+03)   | Discrepancy |
| Sigma2-3  | 1.2e+04 | 2.6e+03 | (9.1e+03 - 1.7e+04) | Discrepancy |
| Sigma2-4  | 1.1e+04 | 2.3e+03 | (7.8e+03 - 1.5e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.2    | Model          |
| CL        | -0.16   | Model          |
| CL        | -0.21   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 4e+03   | Discrepancy    |
| Sigma2-2  | 6.9e+03 | Discrepancy    |
| Sigma2-3  | 1.2e+04 | Discrepancy    |
| Sigma2-4  | 1.1e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-------------------------------------------------
|    |  CL        CL         CL        CL       |
-------------------------------------------------
| CL |  1         -0.057     0.011     0.024    |
| CL |  -0.057    1          -0.065    -0.0078  |
| CL |  0.011     -0.065     1         -0.13    |
| CL |  0.024     -0.0078    -0.13     1        |
-------------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.012      0.045       0.054     |
| Sigma2-2 |  -0.012      1           0.049       -0.05     |
| Sigma2-3 |  0.045       0.049       1           0.012     |
| Sigma2-4 |  0.054       -0.05       0.012       1         |
-------------------------------------------------------------



The MCMC simulation has converged
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
   Duration (HH:MM:SS):                           00:22:40
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.19   | 0.021   | (-0.22 - -0.15)     | Model       |
| CL        | -0.12   | 0.021   | (-0.15 - -0.082)    | Model       |
| CL        | -0.2    | 0.012   | (-0.22 - -0.18)     | Model       |
| CL        | -0.21   | 0.01    | (-0.22 - -0.19)     | Model       |
| Sigma2-1  | 4e+03   | 6.1e+02 | (3.1e+03 - 5.2e+03) | Discrepancy |
| Sigma2-2  | 7.6e+03 | 1.1e+03 | (5.9e+03 - 9.5e+03) | Discrepancy |
| Sigma2-3  | 1.2e+04 | 1.8e+03 | (9.4e+03 - 1.5e+04) | Discrepancy |
| Sigma2-4  | 1.2e+04 | 1.8e+03 | (9.8e+03 - 1.5e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.19   | Model          |
| CL        | -0.12   | Model          |
| CL        | -0.2    | Model          |
| CL        | -0.21   | Model          |
| Sigma2-1  | 4e+03   | Discrepancy    |
| Sigma2-2  | 7.6e+03 | Discrepancy    |
| Sigma2-3  | 1.2e+04 | Discrepancy    |
| Sigma2-4  | 1.2e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
---------------------------------------------
|    |  CL       CL        CL       CL      |
---------------------------------------------
| CL |  1        -0.14     0.06     0.029   |
| CL |  -0.14    1         -0.12    -0.035  |
| CL |  0.06     -0.12     1        -0.2    |
| CL |  0.029    -0.035    -0.2     1       |
---------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.055      -0.018      -0.08     |
| Sigma2-2 |  -0.055      1           0.017       -0.0057   |
| Sigma2-3 |  -0.018      0.017       1           -0.033    |
| Sigma2-4 |  -0.08       -0.0057     -0.033      1         |
-------------------------------------------------------------



The MCMC simulation has converged
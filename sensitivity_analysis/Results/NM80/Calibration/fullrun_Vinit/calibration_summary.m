%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            8751

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            8751

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            8751

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            8751

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
   Duration (HH:MM:SS):                           24:40:18
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.24   | 0.002   | (-0.24 - -0.23)     | Model       |
| CL        | -0.18   | 0.0022  | (-0.18 - -0.17)     | Model       |
| CL        | -0.15   | 0.0071  | (-0.16 - -0.14)     | Model       |
| CL        | -0.21   | 0.00095 | (-0.21 - -0.21)     | Model       |
| Sigma2-1  | 4.1e+03 | 60      | (4e+03 - 4.2e+03)   | Discrepancy |
| Sigma2-2  | 8.1e+03 | 1.2e+02 | (7.9e+03 - 8.3e+03) | Discrepancy |
| Sigma2-3  | 1.5e+04 | 2.3e+02 | (1.4e+04 - 1.5e+04) | Discrepancy |
| Sigma2-4  | 1.2e+04 | 1.8e+02 | (1.2e+04 - 1.2e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Map     | Parameter Type |
----------------------------------------
| CL        | -0.23   | Model          |
| CL        | -0.18   | Model          |
| CL        | -0.15   | Model          |
| CL        | -0.21   | Model          |
| Sigma2-1  | 4.1e+03 | Discrepancy    |
| Sigma2-2  | 8.1e+03 | Discrepancy    |
| Sigma2-3  | 1.5e+04 | Discrepancy    |
| Sigma2-4  | 1.2e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-----------------------------------------------
|    |  CL        CL        CL        CL      |
-----------------------------------------------
| CL |  1         -0.044    -0.013    0.0062  |
| CL |  -0.044    1         -0.18     0.093   |
| CL |  -0.013    -0.18     1         -0.54   |
| CL |  0.0062    0.093     -0.54     1       |
-----------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.037      0.063       -0.008    |
| Sigma2-2 |  -0.037      1           -0.017      0.05      |
| Sigma2-3 |  0.063       -0.017      1           0.022     |
| Sigma2-4 |  -0.008      0.05        0.022       1         |
-------------------------------------------------------------

%------------------- Convergence
   Gelman-Rubin MPSRF:                            1.524059e+00
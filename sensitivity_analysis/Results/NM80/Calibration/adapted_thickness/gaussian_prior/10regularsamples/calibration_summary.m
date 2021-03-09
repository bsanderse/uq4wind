
%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            11

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            11

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            11

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            11

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
   Duration (HH:MM:SS):                           00:11:33
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.11   | 0.095   | (-0.25 - 0.061)     | Model       |
| CL        | -0.094  | 0.076   | (-0.21 - 0.04)      | Model       |
| CL        | -0.22   | 0.053   | (-0.29 - -0.12)     | Model       |
| CL        | -0.2    | 0.043   | (-0.27 - -0.13)     | Model       |
| Sigma2-1  | 1.2e+04 | 6.6e+03 | (4.9e+03 - 2.5e+04) | Discrepancy |
| Sigma2-2  | 1.4e+04 | 7.2e+03 | (6.3e+03 - 2.9e+04) | Discrepancy |
| Sigma2-3  | 2.7e+04 | 1e+04   | (1.2e+04 - 4.6e+04) | Discrepancy |
| Sigma2-4  | 2.3e+04 | 1e+04   | (1e+04 - 4.3e+04)   | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.11   | Model          |
| CL        | -0.094  | Model          |
| CL        | -0.22   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 1.2e+04 | Discrepancy    |
| Sigma2-2  | 1.4e+04 | Discrepancy    |
| Sigma2-3  | 2.7e+04 | Discrepancy    |
| Sigma2-4  | 2.3e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
--------------------------------------------------
|    |  CL         CL        CL         CL       |
--------------------------------------------------
| CL |  1          -0.013    -0.0016    -0.0084  |
| CL |  -0.013     1         0.023      -0.045   |
| CL |  -0.0016    0.023     1          -0.15    |
| CL |  -0.0084    -0.045    -0.15      1        |
--------------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           0.011       -0.015      -0.042    |
| Sigma2-2 |  0.011       1           0.017       -0.0098   |
| Sigma2-3 |  -0.015      0.017       1           0.059     |
| Sigma2-4 |  -0.042      -0.0098     0.059       1         |
-------------------------------------------------------------



The MCMC simulation has converged
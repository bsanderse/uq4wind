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
   Duration (HH:MM:SS):                           00:17:55
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.22   | 0.018   | (-0.25 - -0.19)     | Model       |
| CL        | -0.15   | 0.019   | (-0.18 - -0.12)     | Model       |
| CL        | -0.22   | 0.012   | (-0.24 - -0.2)      | Model       |
| CL        | -0.2    | 0.0091  | (-0.22 - -0.19)     | Model       |
| Sigma2-1  | 3.2e+03 | 4.4e+02 | (2.6e+03 - 4e+03)   | Discrepancy |
| Sigma2-2  | 6.3e+03 | 9.4e+02 | (4.9e+03 - 8e+03)   | Discrepancy |
| Sigma2-3  | 1.2e+04 | 1.7e+03 | (9.8e+03 - 1.5e+04) | Discrepancy |
| Sigma2-4  | 1e+04   | 1.5e+03 | (8e+03 - 1.3e+04)   | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.22   | Model          |
| CL        | -0.15   | Model          |
| CL        | -0.22   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 3.2e+03 | Discrepancy    |
| Sigma2-2  | 6.3e+03 | Discrepancy    |
| Sigma2-3  | 1.2e+04 | Discrepancy    |
| Sigma2-4  | 1e+04   | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
---------------------------------------------
|    |  CL       CL        CL       CL      |
---------------------------------------------
| CL |  1        -0.16     0.014    0.012   |
| CL |  -0.16    1         0.016    0.0083  |
| CL |  0.014    0.016     1        -0.12   |
| CL |  0.012    0.0083    -0.12    1       |
---------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           0.005       -0.013      0.018     |
| Sigma2-2 |  0.005       1           0.041       0.015     |
| Sigma2-3 |  -0.013      0.041       1           0.02      |
| Sigma2-4 |  0.018       0.015       0.02        1         |
-------------------------------------------------------------

%------------------- Convergence
   Gelman-Rubin MPSRF:                            1.533405e+00

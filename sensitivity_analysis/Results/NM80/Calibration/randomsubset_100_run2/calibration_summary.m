
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
   Duration (HH:MM:SS):                           00:16:51
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.22   | 0.016   | (-0.24 - -0.19)     | Model       |
| CL        | -0.16   | 0.019   | (-0.19 - -0.12)     | Model       |
| CL        | -0.14   | 0.059   | (-0.23 - -0.038)    | Model       |
| CL        | -0.2    | 0.0084  | (-0.21 - -0.19)     | Model       |
| Sigma2-1  | 3.2e+03 | 4.7e+02 | (2.5e+03 - 4.1e+03) | Discrepancy |
| Sigma2-2  | 6.4e+03 | 9.4e+02 | (5e+03 - 8.1e+03)   | Discrepancy |
| Sigma2-3  | 1.2e+04 | 1.8e+03 | (9.7e+03 - 1.6e+04) | Discrepancy |
| Sigma2-4  | 1e+04   | 1.4e+03 | (8e+03 - 1.3e+04)   | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.22   | Model          |
| CL        | -0.16   | Model          |
| CL        | -0.14   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 3.2e+03 | Discrepancy    |
| Sigma2-2  | 6.4e+03 | Discrepancy    |
| Sigma2-3  | 1.2e+04 | Discrepancy    |
| Sigma2-4  | 1e+04   | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
---------------------------------------------
|    |  CL        CL       CL       CL      |
---------------------------------------------
| CL |  1         -0.14    0.045    -0.032  |
| CL |  -0.14     1        -0.22    0.098   |
| CL |  0.045     -0.22    1        -0.53   |
| CL |  -0.032    0.098    -0.53    1       |
---------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           0.027       -0.00013    -0.00092  |
| Sigma2-2 |  0.027       1           0.0051      -0.071    |
| Sigma2-3 |  -0.00013    0.0051      1           -0.022    |
| Sigma2-4 |  -0.00092    -0.071      -0.022      1         |
-------------------------------------------------------------



The MCMC simulation has converged

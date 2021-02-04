%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            1000

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            1000

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            1000

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            1000

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
   Duration (HH:MM:SS):                           00:55:20
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.24   | 0.0057  | (-0.25 - -0.23)     | Model       |
| CL        | -0.18   | 0.0067  | (-0.19 - -0.17)     | Model       |
| CL        | -0.13   | 0.022   | (-0.17 - -0.095)    | Model       |
| CL        | -0.2    | 0.0029  | (-0.21 - -0.2)      | Model       |
| Sigma2-1  | 4.1e+03 | 1.8e+02 | (3.8e+03 - 4.4e+03) | Discrepancy |
| Sigma2-2  | 8e+03   | 3.7e+02 | (7.5e+03 - 8.7e+03) | Discrepancy |
| Sigma2-3  | 1.4e+04 | 6.4e+02 | (1.3e+04 - 1.5e+04) | Discrepancy |
| Sigma2-4  | 1.2e+04 | 4.9e+02 | (1.1e+04 - 1.2e+04) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.24   | Model          |
| CL        | -0.18   | Model          |
| CL        | -0.13   | Model          |
| CL        | -0.2    | Model          |
| Sigma2-1  | 4.1e+03 | Discrepancy    |
| Sigma2-2  | 8e+03   | Discrepancy    |
| Sigma2-3  | 1.4e+04 | Discrepancy    |
| Sigma2-4  | 1.2e+04 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-------------------------------------------
|    |  CL       CL       CL       CL     |
-------------------------------------------
| CL |  1        -0.12    0.06     -0.04  |
| CL |  -0.12    1        -0.22    0.17   |
| CL |  0.06     -0.22    1        -0.57  |
| CL |  -0.04    0.17     -0.57    1      |
-------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.018      0.039       0.08      |
| Sigma2-2 |  -0.018      1           -0.023      -0.038    |
| Sigma2-3 |  0.039       -0.023      1           -0.06     |
| Sigma2-4 |  0.08        -0.038      -0.06       1         |
-------------------------------------------------------------


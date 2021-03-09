%----------------------- Inversion output -----------------------%
   Number of calibrated model parameters:         4
   Number of non-calibrated model parameters:     0

   Number of calibrated discrepancy parameters:   4

%------------------- Data and Discrepancy
%  Data-/Discrepancy group 1:
   Number of independent observations:            50

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       1

%  Data-/Discrepancy group 2:
   Number of independent observations:            50

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            50

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            50

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
   Duration (HH:MM:SS):                           00:16:43
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.23   | 0.017   | (-0.26 - -0.2)      | Model       |
| CL        | -0.17   | 0.014   | (-0.19 - -0.15)     | Model       |
| CL        | -0.23   | 0.0087  | (-0.24 - -0.21)     | Model       |
| CL        | -0.21   | 0.007   | (-0.22 - -0.2)      | Model       |
| Sigma2-1  | 1.5e+03 | 3.1e+02 | (1.1e+03 - 2.1e+03) | Discrepancy |
| Sigma2-2  | 2e+03   | 4.3e+02 | (1.4e+03 - 2.8e+03) | Discrepancy |
| Sigma2-3  | 3.6e+03 | 8e+02   | (2.4e+03 - 5e+03)   | Discrepancy |
| Sigma2-4  | 3e+03   | 6.7e+02 | (2.1e+03 - 4.2e+03) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.23   | Model          |
| CL        | -0.17   | Model          |
| CL        | -0.23   | Model          |
| CL        | -0.21   | Model          |
| Sigma2-1  | 1.5e+03 | Discrepancy    |
| Sigma2-2  | 2e+03   | Discrepancy    |
| Sigma2-3  | 3.6e+03 | Discrepancy    |
| Sigma2-4  | 3e+03   | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
----------------------------------------------
|    |  CL        CL        CL       CL      |
----------------------------------------------
| CL |  1         -0.034    0.028    -0.041  |
| CL |  -0.034    1         -0.03    -0.026  |
| CL |  0.028     -0.03     1        -0.16   |
| CL |  -0.041    -0.026    -0.16    1       |
----------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           0.019       0.072       -0.0092   |
| Sigma2-2 |  0.019       1           0.05        -0.055    |
| Sigma2-3 |  0.072       0.05        1           -0.044    |
| Sigma2-4 |  -0.0092     -0.055      -0.044      1         |
-------------------------------------------------------------




R_hat_full =

    1.5404
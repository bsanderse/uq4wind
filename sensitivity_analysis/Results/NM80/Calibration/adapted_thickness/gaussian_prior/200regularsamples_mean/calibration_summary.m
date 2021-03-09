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
   Number of independent observations:            200

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       2

%  Data-/Discrepancy group 3:
   Number of independent observations:            200

   Discrepancy:
      Type:                                       Gaussian
      Discrepancy family:                         Scalar
      Discrepancy parameters known:               No

   Associated outputs:
      Model 1: 
         Output dimensions:                       3

%  Data-/Discrepancy group 4:
   Number of independent observations:            200

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
   Duration (HH:MM:SS):                           00:24:12
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.23   | 0.011   | (-0.25 - -0.22)     | Model       |
| CL        | -0.17   | 0.011   | (-0.19 - -0.15)     | Model       |
| CL        | -0.23   | 0.0071  | (-0.24 - -0.21)     | Model       |
| CL        | -0.21   | 0.0056  | (-0.22 - -0.2)      | Model       |
| Sigma2-1  | 2.9e+03 | 2.9e+02 | (2.4e+03 - 3.4e+03) | Discrepancy |
| Sigma2-2  | 4.8e+03 | 5.1e+02 | (4.1e+03 - 5.7e+03) | Discrepancy |
| Sigma2-3  | 9.1e+03 | 9.4e+02 | (7.7e+03 - 1.1e+04) | Discrepancy |
| Sigma2-4  | 7.2e+03 | 7.3e+02 | (6e+03 - 8.5e+03)   | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.23   | Model          |
| CL        | -0.17   | Model          |
| CL        | -0.23   | Model          |
| CL        | -0.21   | Model          |
| Sigma2-1  | 2.9e+03 | Discrepancy    |
| Sigma2-2  | 4.8e+03 | Discrepancy    |
| Sigma2-3  | 9.1e+03 | Discrepancy    |
| Sigma2-4  | 7.2e+03 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-----------------------------------------------
|    |  CL        CL        CL        CL      |
-----------------------------------------------
| CL |  1         -0.05     -0.068    0.067   |
| CL |  -0.05     1         0.018     -0.011  |
| CL |  -0.068    0.018     1         -0.18   |
| CL |  0.067     -0.011    -0.18     1       |
-----------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.076      -0.013      0.031     |
| Sigma2-2 |  -0.076      1           0.024       0.0029    |
| Sigma2-3 |  -0.013      0.024       1           0.032     |
| Sigma2-4 |  0.031       0.0029      0.032       1         |
-------------------------------------------------------------




R_hat_full =

    1.5034
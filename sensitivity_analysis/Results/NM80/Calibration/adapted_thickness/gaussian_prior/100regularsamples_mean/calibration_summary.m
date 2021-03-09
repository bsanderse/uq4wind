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
   Duration (HH:MM:SS):                           00:24:14
   Number of sample points:                       1.00e+05

%------------------- Posterior Marginals
---------------------------------------------------------------------
| Parameter | Mean    | Std     | (0.05-0.95) Quant.  | Type        |
---------------------------------------------------------------------
| CL        | -0.23   | 0.013   | (-0.25 - -0.21)     | Model       |
| CL        | -0.17   | 0.013   | (-0.19 - -0.15)     | Model       |
| CL        | -0.22   | 0.0077  | (-0.24 - -0.21)     | Model       |
| CL        | -0.21   | 0.0064  | (-0.22 - -0.2)      | Model       |
| Sigma2-1  | 2e+03   | 3e+02   | (1.6e+03 - 2.5e+03) | Discrepancy |
| Sigma2-2  | 3.2e+03 | 4.6e+02 | (2.5e+03 - 4e+03)   | Discrepancy |
| Sigma2-3  | 6e+03   | 8.6e+02 | (4.7e+03 - 7.5e+03) | Discrepancy |
| Sigma2-4  | 4.7e+03 | 6.8e+02 | (3.7e+03 - 5.9e+03) | Discrepancy |
---------------------------------------------------------------------

%------------------- Point estimate
----------------------------------------
| Parameter | Mean    | Parameter Type |
----------------------------------------
| CL        | -0.23   | Model          |
| CL        | -0.17   | Model          |
| CL        | -0.22   | Model          |
| CL        | -0.21   | Model          |
| Sigma2-1  | 2e+03   | Discrepancy    |
| Sigma2-2  | 3.2e+03 | Discrepancy    |
| Sigma2-3  | 6e+03   | Discrepancy    |
| Sigma2-4  | 4.7e+03 | Discrepancy    |
----------------------------------------

%------------------- Correlation matrix (Model Parameters)
-----------------------------------------------
|    |  CL        CL        CL        CL      |
-----------------------------------------------
| CL |  1         -0.064    0.032     0.065   |
| CL |  -0.064    1         -0.038    0.0026  |
| CL |  0.032     -0.038    1         -0.15   |
| CL |  0.065     0.0026    -0.15     1       |
-----------------------------------------------

%------------------- Correlation matrix (Discrepancy Parameters)
-------------------------------------------------------------
|          |  Sigma2-1    Sigma2-2    Sigma2-3    Sigma2-4  |
-------------------------------------------------------------
| Sigma2-1 |  1           -0.0037     -0.049      0.0081    |
| Sigma2-2 |  -0.0037     1           -0.039      0.032     |
| Sigma2-3 |  -0.049      -0.039      1           0.023     |
| Sigma2-4 |  0.0081      0.032       0.023       1         |
-------------------------------------------------------------




R_hat_full =

    2.577
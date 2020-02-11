function Y = almanac_vectorized(Q,A)
% B: uncertain parameter, size 1*p
% A: design matrix, contains independent variable, size N*p: 

% UQLab uses convention that y is a row vector if a single beta (row vector) is given
% when vectorized is true, it can then determine a full matrix of Y for a
% given matrix of beta
% this means that instead of A*B we write B*A'

Y = Q*A';
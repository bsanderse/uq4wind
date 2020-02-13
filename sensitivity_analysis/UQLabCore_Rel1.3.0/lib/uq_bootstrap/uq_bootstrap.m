function resIDX = uq_bootstrap(N,B)

%% consistency check: X must be a 2D matrix 
resIDX = round(rand(B,N)*(N-1)) + 1;
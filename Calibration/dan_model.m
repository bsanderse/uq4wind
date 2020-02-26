function Y = dan_model(Q,beta)
% Q: uncertain parameter, size p*2
% beta: contains independent variable, size N_data*p: 

% y = b_1 + b_2*x 
Y = Q(1).*sin(beta) + Q(2).*cos(beta);


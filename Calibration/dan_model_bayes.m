function Y = dan_model_bayes(X)
% UQ_SIMPLYSUPPORTEDBEAM computes the midspan deflection of a simply 
% supported beam under uniform loading
%
%   Model with five input parameters  X= [beta c W cL cD]
%      beta:   Flow angle
%         c:   chord
%         W:   Relative velocity      
%        cL:   Lift coefficient
%        cD:   Drag coefficient
%         
%
%   Output:  Y = 0.5*1.225*(W^2)*c*(cL*sin(beta) + cD*sin(beta))

% Vectorized implementation
Y = 0.5*1.225.*X(:, 2).*(X(:, 3).^2).*(X(:, 4).*sin(X(:, 1)) + X(:, 5).*cos(X(:, 1)));
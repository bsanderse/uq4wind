function Y = linearmodel_vectoroutput(B)
% B: uncertain parameter, size 1*p, or M*p
% p: dimension of parameter vector
% M: number of parameter vectors to be evaluated

% output Y: size M 

Y = B*[1 3; -5 6]';


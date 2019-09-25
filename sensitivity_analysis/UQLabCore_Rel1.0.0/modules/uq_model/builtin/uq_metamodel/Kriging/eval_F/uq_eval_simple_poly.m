function VALUE = uq_eval_simple_poly(ORDER, X)
% VALUE = UQ_EVAL_SIMPLE_POLY(ORDER, X): evaluate the univariate
%     polynomials of the form 1 + X + X^2 + ... + X^ORDER on the vector X.

    
% make sure X is a column vector, and not a row vector
if size(X,2) ~= 1
   error('uq_eval_simple_poly is designed to work with X in column vector format');
end


%% RECURSIVELY CALCULATE THE POLYNOMIAL VALUE UP TO THE REQUESTED ORDER
VALUE = ones(numel(X), ORDER + 1);
for ii = 1 : ORDER
   VALUE(:,ii+1) = VALUE(:,ii).*X;
end





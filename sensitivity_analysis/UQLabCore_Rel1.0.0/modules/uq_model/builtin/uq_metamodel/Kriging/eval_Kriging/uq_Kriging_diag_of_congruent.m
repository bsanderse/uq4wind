function Diag = uq_Kriging_diag_of_congruent( A,B )
% Diag = UQ_KRIGING_DIAG_OF_CONGRUENT(A,B): Returns the diagonal elements of a congruent
% matrix of a special form.
%   Diag = diag(C), where C = A' * B^(-1) * A

if size(A,1) ~= size(B,2)
    A = transpose(A) ;
end
% For returning a row vector that contains the diagonal elements of C:
if rcond(B) > eps
    % for better numerical precision
     C = B \ A ;  
     Diag = sum(A.*C, 1);
else
    Binv = pinv(B) ;
    Diag = sum(( Binv * A).* A , 1) ;
end

function xAug = uq_inversion_AugmentInput(x,nonConst,idConst,valConst)
% UQ_INVERSION_AUGMENTINPUT augments a sample with constant values
%   
%   XAUG = UQ_INVERSION_AUGMENTINPUT(X, NONCONST, IDCONST, VALCONST)
%   augments the sample in X with the constant values in VALCONST
%
%   Additional notes:
%
%   - This function is only required if the prior distribution contains
%     constants and then serves as a wrapper before the forward model is
%     evaluated.

% Initialize
[m,n] = size(x);
nAug = n + length(idConst);

% pre-allocation
xAug = zeros(m,nAug);

% loop over rows
j = 0;							% running index for variables
k = 0;							% running index for constants
for i = 1:nAug
    if any(i==nonConst)						% variable parameters
      j = j + 1;
      xAug(:,i) = x(:,j);
    elseif any(i==idConst)					% constant parameters
      k = k + 1;
      xAug(:,i) = ones(m,1) * valConst(k);
    end
end
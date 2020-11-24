function [maxnp, idx] = maxnonpositive(Array, varargin)

switch nargin
    
    case 1
        
    Array(Array>0) = nan;
    
    [maxnp, idx] = max(Array);
    
    otherwise
        
        disp('Invalid number of arguments')
        
end
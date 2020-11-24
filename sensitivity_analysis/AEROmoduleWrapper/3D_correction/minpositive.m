function [minp, idx] = minpositive(Array, varargin)

switch nargin
    
    case 1
        
    Array(Array<=0) = nan;
  
    [minp, idx] = min(Array);
    
    otherwise
        
        disp('Invalid number of arguments')
        
end
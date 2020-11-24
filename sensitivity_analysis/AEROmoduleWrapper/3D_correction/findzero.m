function [fzero] = findzero(x, f, varargin)

%estimate the value of a function at 0, from given (Arrays) function values at
%different points

switch nargin 
    
    case 2
        
    [maxnp, idx1] = maxnonpositive(x);
    
    if isnan(maxnp)
        
        disp('No f values when x<=0, cannot interpolate')
        return
        
    end
    
    [minp, idx2] = minpositive(x);
    
    if isnan(minp)
        
        disp('No f values when x>0, cannot interpolate')
        return
        
    end
    
    if maxnp == 0
        
        idx = idx1;
        fzero = f(idx);
        
    else
        
        xval = [maxnp minp];
        
        fval = [f(idx1) f(idx2)];
        
        mval = [maxnp 0 minp];
        
        full_f = interp1(xval, fval, mval,'linear');
        
        fzero = full_f(2);
        
    end
    
    otherwise
        
        disp('Invalid number of arguments')
end
    

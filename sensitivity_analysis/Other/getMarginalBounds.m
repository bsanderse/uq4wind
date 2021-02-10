function [ domain ] = getMarginalBounds( Marginal )
%GETMARGINALBOUNDS get suitable bounds to display a PDF
%   Marginal is an object following the UQLab standard

switch Marginal.Type
    case 'Gaussian'
        domain = [Marginal.Moments(1)-2*Marginal.Moments(2), ...
                  Marginal.Moments(1)+2*Marginal.Moments(2)];
    case 'Uniform'
        domain = [Marginal.Bounds(1),Marginal.Bounds(2)];
    case 'Weibull'
        domain = [0,Marginal.Moments(1) + 4*Marginal.Moments(2)];
    otherwise
        error('unknown pdf type in getMarginalBounds');
end

end


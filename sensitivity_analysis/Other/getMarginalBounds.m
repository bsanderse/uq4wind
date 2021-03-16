function [ domain ] = getMarginalBounds( Marginal )
%GETMARGINALBOUNDS get suitable bounds to display a PDF
%   Marginal is an object following the UQLab standard

switch Marginal.Type
    case 'Gaussian'
        if (~isfield(Marginal,'Bounds'))
            domain = [Marginal.Moments(1)-3*Marginal.Moments(2), ...
                      Marginal.Moments(1)+3*Marginal.Moments(2)];
        else
            domain = [Marginal.Bounds(1),Marginal.Bounds(2)];
        end
    case 'Uniform'
        if (~isfield(Marginal,'Bounds'))
            domain = [Marginal.Parameters(1),Marginal.Parameters(2)];
        else
            domain = [Marginal.Bounds(1),Marginal.Bounds(2)];
        end

    case 'Weibull'
        domain = [0,Marginal.Moments(1) + 4*Marginal.Moments(2)];
    otherwise
        error('unknown pdf type in getMarginalBounds');
end

end


function [pdf,x] = truncnormpdf(x,mu,sigma,a,b)
% pdf of truncated normal distribution with parameters mu and sigma,
% truncated to the interval a, b

% non-truncated pdf:
num  = normpdf(x,mu,sigma);
% truncate by scaling with the area as given by the cdf:
den2 = normcdf(b,mu,sigma);
den1 = normcdf(a,mu,sigma);

pdf  = (x>a & x<b).* num /(den2 - den1);


end





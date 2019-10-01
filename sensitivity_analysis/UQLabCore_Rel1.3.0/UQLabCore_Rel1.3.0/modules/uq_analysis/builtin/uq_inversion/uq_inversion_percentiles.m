function percentiles = uq_inversion_percentiles(Sample3D, probabilities)
% UQ_INVERSION_PERCENTILES computes the percentiles specified by
%     PROBABILITIES using the passed SAMPLE3D. 
%
% See also UQ_POSTPROCESSINVERSION, UQ_INVERSION_MEAN

% reshape Sample3D to Sample2D
nDim = size(Sample3D,2);
Sample2D = reshape(permute(Sample3D,[2 1 3]),nDim,[])';

% compute percentiles
percentiles = quantile(Sample2D,probabilities);
end


function [Fcomplex,Fabs,Fphi] = getFourierCoefficientsRealData(F,n_fourier)
% F is supposed to be a real matrix of size n*n_col,
% and an FFT is performed for each column of F
% the FFT coefficients are sorted according to their PSD
% the mean is always reported as first entry, independent whether is has a
% large PSD

n     = size(F,1);
n_col = size(F,2);

% n_fourier is the number of Fourier modes to keep
if (nargin == 1)
    n_fourier = floor(n/2); % report all coefficients
end

% actual number of coefficients is 2*n_fourier - 1 
Fcomplex = zeros(2*n_fourier-1,n_col);

for k = 1:n_col
    
        % do the Fourier transform
    % include scaling to get physical interpretable coefficients
    Fhat = fft(F(:,k),n)/n;
    
    % the first entry is the constant mode (mean)
    % this is not necessarily the maximum power spectral density
    % but we want this in any case to be part of the output
    ind_mean  = 1;
    Fhat_mean = Fhat(ind_mean);
    
    % remove mean
    Fhat(ind_mean) = [];
    
    % get the power spectral density of the remaining terms
    PSD  = Fhat.*conj(Fhat);
    
    % select indices with largest PSD by sorting the PSD
    [~,ind] = sort(PSD,'desc');
    
    %
    Fcomplex(:,k) = [Fhat_mean;Fhat(ind(1:2*n_fourier-2))];
   
end

Fabs = abs(Fcomplex);
Fphi = angle(Fcomplex);

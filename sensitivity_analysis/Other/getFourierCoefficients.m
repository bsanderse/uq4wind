function [Fcomplex,Fabs,Fphi] = getFourierCoefficients(F,n_fourier)
% F is supposed to be a real matrix of size n*n_col,
% and an FFT is performed for each column of F

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
    
    % get the power spectral density
    PSD  = Fhat.*conj(Fhat);

    % first half of frequencies contains all information, because the signal is
    % real, so c_k = c_{-k}
    % we therefore plot the one-sided (positive) frequency range only
%     freqVals = (1/dt)*(0:floor(n/2)-1)'/n;
    
    % power in time and in frequency domain (should match)
%     norm(Fn_k)^2/n;
%     sum(PSD);

    % for plotting purposes, the fftshift can be useful:
    % Fhat_shift = fftshift(Fhat);
    % PSD_plot = Fhat_shift.*conj(Fhat_shift);
    % full frequency range:
    % freq = 1/(dt*n)*(0:n);
    % freq = (1/dt)*(-n/2:n/2-1)/n;
    
    % select indices with largest PSD by sorting the PSD
    [~,ind] = sort(PSD,'desc');
    %
    ind_select = ones(n,1);
    % select 2*n_keep indices, where the factor 2 is needed because we need the
    % coefficients associated with negative frequencies as well to do the inverse
    % transform
%     ind_select(ind(2*n_fourier:end)) = 0;
    % set other indices to 0
%     Fout(:,k) = ind_select.*Fhat;
    % back to the time domain
%     Fnew     = n*ifft(Fhat_new,n);
    
    % note: if we want physical meaning out of the Fourier coefficients
    % (amplitude, angle) we need to multiply by 2 if we don't
    % include the complex conjugates
    % below is with all coefficients:
    ind_keep  = ind(1:2*n_fourier-1);
    Fcomplex(:,k) = Fhat(ind_keep);
    

%     abs(f_keep);
%     angle(f_keep);

    % alternatively, with only positive frequencies:
%     ind_pos = ind(2:2:2*(n_keep-1));
%     f_pos  = Fhat(ind_pos);
%     abs(Fhat(ind(1)))
%     abs(f_pos)*2
%     angle(f_pos)    
    
    
%     figure(9)
%     semilogy(freqVals(1:end),PSD(1:floor(n/2)),'x-','Color',colormap(k,:))
%     hold on
%     semilogy(freqVals.*ind_select(1:floor(n/2)),PSD(1:floor(n/2)).*ind_select(1:floor(n/2)),'s','Color',colormap(k,:))
%     
%     figure(10)
%     plot(azi_exp_data,Fn_pert,'-','Color',colormap(k,:));
%     hold on
%     plot(azi_exp_data,Fnew,'--','Color',colormap(k,:));
    
end

Fabs = abs(Fcomplex);
Fphi = angle(Fcomplex);

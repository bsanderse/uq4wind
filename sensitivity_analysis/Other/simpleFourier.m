%% simple example of  FFT and iFFT
% https://www.gaussianwaves.com/2014/07/how-to-plot-fft-using-matlab-fft-of-basic-signals-sine-and-cosine-waves/
% https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/

n  = 100; % number of samples
t  = linspace(0,1,n)';
dt = t(2)-t(1);
freq = 1/dt; %sampling frequency should be at least 2*highest frequency in signal
f    = 2 + 0.1*cos(2*pi*5*t + 10*pi/180);% + sin(2*pi*120*t);
n    = length(t);
fhat = fft(f,n)/n; % if scaling is included it needs to be added also in ifft
psd  = fhat.*conj(fhat); %/(n^2); % n^2 disappears if it we do fft/n

% power in time and in frequency domain (should match)
norm(f)^2/n;
sum(psd);

fVals = freq*(0:n/2-1)/n;
figure
semilogy(fVals,psd(1:n/2))

% select indices with largest PSD by sorting the PSD
[val,ind] = sort(psd,'desc');
%
ind_select = ones(n,1);
% select 2*n_keep indices, where the factor 2 is needed because we need the
% coefficients associated with negative frequencies as well to do the inverse
% transform
n_keep = 3;
ind_select(ind(2*n_keep:end)) = 0;
% set other indices to 0
fhat_new = ind_select.*fhat;
fnew    = n*ifft(fhat_new,n);

% note: if we want physical meaning out of the Fourier coefficients
% (amplitude, angle) we need to divide by n (and multiply by 2 if we don't
% include the complex conjugates)
ind_keep = ind(1:2*n_keep-1);
f_keep = fhat(ind_keep);
abs(f_keep);

% the phase angle is more tricky business
% https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/
%angle(f_keep)
% f2 = f_keep;
% threshold = max(abs(f_keep))/10000; %tolerance threshold
% f2(abs(f_keep)<threshold) = 0; %maskout values that are below the threshold
% phase=atan2(imag(f2),real(f2))*180/pi; %phase information

% get only positive frequencies and mean
ind_pos = 2:2:2*(n_keep-1);
abs(fhat(1))
abs(fhat(ind_keep(ind_pos)))*2
f2 = fhat;
threshold = max(abs(f2))/1e4; %tolerance threshold
f2(abs(f2)<threshold) = 0; %maskout values that are below the threshold
angle(f2(ind_keep(ind_pos)))*180/pi


figure
plot(t,f);
hold on
plot(t,fnew);
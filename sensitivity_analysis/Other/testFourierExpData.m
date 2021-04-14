%% test Fourier representation of azimuth-dependent data
clearvars
% close all


%% figure numbers
figID_PSD = 9;
figID_Fourier = 105;

%% load  experimental data
% run this file from sensitivity_analysis folder
RPM = 425.1;

root_folder = pwd;

% NewMexicoData as obtained from Koen Boorsma (TNO)
folder_exp    = fullfile(root_folder,'..','Experimental','NewMexicoData');
filename_exp = 'R52P83D948_loads_n.dat';
full_filename_exp = fullfile(folder_exp,filename_exp);
% read in the table
output_raw   = readNewMexicoModified(full_filename_exp);    

% the position of the sections of the experimental data which are used for
% interpolation of the aeromodule results: see NM80_calibrate_readoutput.m
r_exp_data = [0.25 0.35 0.6 0.82 0.92]*2.25;

azi_exp_data = output_raw.Azi;

% Because the model has different discrepancy options at different radial locations,
% the measurement data is stored in five different data structures:
% normal forces at five stations
Fn_exp_data = table2array(output_raw(:,4:8));

% interpolate the data to be equidistant in time/azimuth
dazi  = 10;
dt    = dazi/RPM/6;
azi_exp_int = (0:dazi:360-dazi)';
Fn_exp_int  = spline(azi_exp_data',Fn_exp_data',azi_exp_int)';



% make plot of data
% figure(10)
% hold on
% plot(azi_exp_data,Fn_exp_data,'x-');
% plot(azi_exp_int,Fn_exp_int,'s-');

%% loop over radial sections and do FFT for each section

figure(figID_Fourier)
colormap = get(gca,'ColorOrder');

% radial indices to plot:
r_index = 1:5;

% number of fourier modes to keep (including the mean)
% n_keep = 1 :  zeroth mode (meanonly) (1 coefficient in total)
% n_keep = 2 :  + amplitude and phase of first mode (3 coefficients in
% total)
% n_keep = 3 :  + amplitude and phase of second mode (5 coefficients in
% total)
n_keep = 2;

%     

for k = 1:length(r_index)
    
    
    %% Fourier
    n    = length(azi_exp_data);
%     dazi = mean(diff(azi_last_rev)); % in degrees
%     dt   = mean(diff(t_last_rev)); % in seconds
    Fn_k   = Fn_exp_data(:,r_index(k));
    % subtract the mean of the data (possible, not required)
%     Fn_mean = mean(Fn);
    Fn_pert = Fn_k; % - Fn_mean;
    % do the fourier transform
    % include scaling to get physical interpretable coefficients
    Fhat = fft(Fn_pert,n)/n;
    % get the power spectral density
    PSD  = Fhat.*conj(Fhat);

    % first half of frequencies contains all information, because the signal is
    % real, so c_k = c_{-k}
    % we therefore plot the one-sided (positive) frequency range only
    freqVals = (1/dt)*(0:floor(n/2)-1)'/n;
    
    % power in time and in frequency domain (should match)
    norm(Fn_k)^2/n;
    sum(PSD);

    % for plotting purposes, the fftshift can be useful:
    % Fhat_shift = fftshift(Fhat);
    % PSD_plot = Fhat_shift.*conj(Fhat_shift);
    % full frequency range:
    % freq = 1/(dt*n)*(0:n);
    % freq = (1/dt)*(-n/2:n/2-1)/n;
    
    % select indices with largest PSD by sorting the PSD
    [val,ind] = sort(PSD,'desc');
    %
    ind_select = ones(n,1);
    % select 2*n_keep indices, where the factor 2 is needed because we need the
    % coefficients associated with negative frequencies as well to do the inverse
    % transform
    ind_select(ind(2*n_keep:end)) = 0;
    % set other indices to 0
    Fhat_new = ind_select.*Fhat;
    % back to the time domain
    Fnew     = n*ifft(Fhat_new,n);
    
    % note: if we want physical meaning out of the Fourier coefficients
    % (amplitude, angle) we need to multiply by 2 if we don't
    % include the complex conjugates
    % below is with all coefficients:
    ind_keep = ind(1:2*n_keep-1);
    f_keep = Fhat(ind_keep);
    abs(f_keep);
    angle(f_keep);

    % alternatively, with only positive frequencies:
    ind_pos = ind(2:2:2*(n_keep-1));
    f_pos  = Fhat(ind_pos);
    abs(Fhat(ind(1)));
    abs(f_pos)*2;
    angle(f_pos);   
    
    
    figure(figID_PSD)
    semilogy(freqVals(1:end),PSD(1:floor(n/2)),'x-','Color',colormap(k,:))
    hold on
    semilogy(freqVals.*ind_select(1:floor(n/2)),PSD(1:floor(n/2)).*ind_select(1:floor(n/2)),'s','Color',colormap(k,:))
    
    figure(figID_Fourier)
    plot(azi_exp_data,Fn_pert,'s-','Color',colormap(k,:), 'LineWidth',1);
    hold on
%     plot(azi_exp_data,Fnew,'--','Color',colormap(k,:), 'LineWidth',2);
    
end

figure(figID_PSD)
grid on
legend('Section 1','Section 1 - selected modes',...
    'Section 2','Section 2 - selected modes',...
    'Section 3','Section 3 - selected modes',...
    'Section 4','Section 4 - selected modes',...
    'Section 5','Section 5 - selected modes');
xlabel('Frequency [1/s]');
ylabel('Power spectral density');

figure(figID_Fourier)
grid on
legend('Section 1','Section 1 - Fourier',...
    'Section 2','Section 2 - Fourier',...
    'Section 3','Section 3 - Fourier',...
    'Section 4','Section 4 - Fourier',...
    'Section 5','Section 5 - Fourier');
xlabel('azimuth [degree]')
ylabel('F_N [N/m]');
xlim([0 360])
ylim([0 400])
% title('Approximation of experimental data with Fourier modes')




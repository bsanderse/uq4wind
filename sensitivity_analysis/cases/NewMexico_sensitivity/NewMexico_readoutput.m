function Y = NewMexico_readoutput(output_dir,P)

switch P.FixedParameters.QoI

    % return output depending on the QoI we are interested in
    case 'Power'
        filename = fullfile(output_dir,'AeroPower.dat'); % location of the AeroPower.dat output file
        [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename); %, P{23}, P{24});
        Y = mean(PowerWatt);
        
    case 'Axial_Force' % axial force of the rotor
        filename = fullfile(output_dir,'AeroPower.dat'); % location of the AeroPower.dat output file
        [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename); %, P{23}, P{24});
        Y = mean(Axial_ForceN);
        
        
        % (trial)
    case 'Axial_Force_Blade'
        % in this case, the cumulative axial force of a blade (1,2 or 3) will be read 
        filename = fullfile(output_dir,'AeroOutput.txt'); % location of the AeroPower.dat output file
        V = readtable(filename,"ReadVariableNames",true,...
            "PreserveVariableNames",true); %, P{23}, P{24});
        f = V{:,12};
        f_mean = mean(f);
        Y = f_mean;
        
        
    case 'Sectional_normal_force'   
        % in this case the QoI is a vector, returning the time-averaged force at
        % each section
        filename = fullfile(output_dir,'B1n_BEM.txt');
        
        % Read output data generated by AeroModule
        D = readtable(filename,'HeaderLines', 4,"ReadVariableNames",true,...
            "PreserveVariableNames",true); % Reads the variable data from the 
                                           % specified .txt file
        % column 1 is time, 
        % column 2 is azimuth,
        % columns 3:end correspond to different radial locations
        Fn = D{:,3:end};     
        
        % Radial stations   
        r_sim = str2double(D.Properties.VariableNames(3:end)); 
        % the radial stations are expressed in % of the blade length,
        % excluding the hub radius
        % to compare to experiments, add the hub radius:
        r_sim = r_sim + 0.21;
        
        % Locations of experimental data
        % These values are made available from NewMexico: 
        r_exp = P.FixedParameters.r_exp;
        
        switch P.FixedParameters.QoI_type

            case 'mean'
                Fn_mean = mean(Fn,1); % Mean (average in time) values at different radial stations
                
                % Interpolation
                % These values are made available from NewMexico: 
                Y   = spline(r_sim,Fn_mean,r_exp); % Interpolated data using spline

            case 'full'
                
                % number of revolutions to consider (counting backward)
                n_rev     = P.FixedParameters.n_rev;
                % number of Fourier coefficients to keep (including mean)
                % note: we get (n_fourier-1)*2 + 1 coefficients
%                 n_fourier = P.FixedParameters.n_fourier;
                index_fourier = P.FixedParameters.index_fourier;
                fourier_type  = P.FixedParameters.fourier_type;
                % radial indices to consider:
                r_index   = P.FixedParameters.r_index;
                                               
                % Use full (azimuth dependent) solution
                % select a couple of revolutions by looking at where azimuth is smaller
                % than a threshold value; the threshold value is taken as the minimum of the difference
                % between azimuth angles
                t_sim         = D{:,1};
                azi_sim       = D{:,2};
                delta_azi     = floor(min(abs(diff(azi_sim)))); % this should be around 10;
                ind_small_azi = find(azi_sim<delta_azi,n_rev+1,'last');
                % index of requested revolutions
                ind_last_rev  = ind_small_azi(1):ind_small_azi(end)-1;
                
                % select force, azimuth and time based on this index
                Fn_last_rev   = Fn(ind_last_rev,:);
                azi_last_rev  = azi_sim(ind_last_rev);

                % get the time step from the simulation data, it should
                % be constant
                t_last_rev    = t_sim(ind_last_rev);
                dt   = mean(diff(t_last_rev)); % in seconds

                % number of data points
                % n    = length(azi_last_rev); 
                
                % Interpolation: columns of Fn_last_rev are interpolated to
                % yield new columns at r_exp positions
                Fn_int   = spline(r_sim,Fn_last_rev,r_exp);        
                
                % get the coefficients of the dominant modes
                % the coefficients are ordered according to the PSD                
                Fhat     = getFourierCoefficientsRealData(Fn_int); %,n_fourier);
%                 ind_select = 2:2:2*(n_fourier-1);

                
                Y = [];
                n_r_index = length(r_index);
                % loop over radial indices
                for k=1:n_r_index
                    
                    % loop over modes
                    % the convention is that index=1 indicates the mean
                    % index=2 indicates amplitude first mode
                    % index=3 indicates phase shift first mode
                    % index=4 indicates amplitude second mode
                    % etc.
                    for mode = 1:length(index_fourier)
                        ind_select = index_fourier(mode);
                        Fcurr      = Fhat(ind_select,r_index(k));
                        
                        % exception for the mean:
                        if (ind_select==1) % this means the mean is requested
                            F_add = abs(Fcurr);
                        else
                            switch fourier_type
                                case 'amp_phase'
                                    % even index => amplitude
                                    if (mod(ind_select,2)==0)
                                        F_add = 2*abs(Fcurr);
                                    else % odd index => angle
                                        F_add = angle(Fcurr);
                                    end
                                case 'real_imag'
                                    % even index => real part
                                    if (mod(ind_select,2)==0)
                                        F_add = 2*real(Fcurr);
                                    else % odd index => imaginary part
                                        F_add = 2*imag(Fcurr);
                                    end
                                otherwise
                                    error('wrong specification of Fourier type');
                            end
                        end
                        
                        Y = [Y F_add];
                    end
                    
                    % add mean separately
%                     Fhat_mean = abs(Fhat(1,r_index(k)));
%                     Fhat_new  = Fhat(ind_select,r_index(k));
                    
                    % save the complex coefficients in terms of amplitude
                    % and phase angle
                    % since we only store the positive frequencies, we need
                    % to multiply by 2 for the physical amplitudes
%                     Y = horzcat(Y,[Fhat_mean 2*abs(Fhat_new)' angle(Fhat_new)']);
                    
                end
                
%                 figure(100)
%                 plot(t_last_rev,Fn_int(:,1),'-');
%                 hold on
%                 figure(101)
%                 plot(t_last_rev,Fn_int(:,2),'-');
%                 hold on
%                 figure(102)
%                 plot(t_last_rev,Fn_int(:,3),'-');
%                 hold on
%                 figure(103)
%                 plot(t_last_rev,Fn_int(:,4),'-');
%                 hold on
%                 figure(104)
%                 plot(t_last_rev,Fn_int(:,5),'-');
%                 hold on                
                
                
                    
%                 for k=1:length(r_index)
% 
% 
%                     Fn_k = Fn_int(:,r_index(k));
%                     % do the fourier transform
%                     Fhat = fft(Fn_k,n)/n; % include scaling with 1/n to get physical results for coefficients
%                     % get the power spectral density
%                     PSD  = Fhat.*conj(Fhat);
% 
%                     % first half of frequencies contains all information, because the signal is
%                     % real, so c_k = c_{-k}
%                     % we therefore plot the one-sided (positive) frequency range only
% %                     freqVals = (1/dt)*(0:floor(n/2)-1)'/n;
% 
%                     % select indices with largest PSD by sorting the PSD
%                     [~,ind] = sort(PSD,'desc');
%                     %
% %                     ind_select = ones(n,1);
%                     % select 2*n_keep indices, where the factor 2 is needed because we need the
%                     % coefficients associated with negative frequencies as well to do the inverse
%                     % transform
% %                     ind_select(ind(2*n_keep:end)) = 0; 
%                     
%                     % alternatively, we can select directly from Fhat:
%                     % note that the indices that are skipped correspond to
%                     % the complex conjugate, so they don't need to be
%                     % stored
%                     ind_select = 2:2:2*(n_fourier-1);
%                     
%                     % add mean separately
%                     Fhat_mean = abs(Fhat(ind(1)));
%                     Fhat_new  = Fhat(ind(ind_select));
% 
%                     % save the complex coefficients in terms of amplitude
%                     % and phase angleei
%                     % since we only store the positive frequencies, we need
%                     % to multiply by 2 for the physical amplitudes
%                     Y = horzcat(Y,[Fhat_mean 2*abs(Fhat_new)' angle(Fhat_new)']);
%                     
%                 end
        
        end
        
    otherwise
        error(strcat('QoI type unknown'));
        
end
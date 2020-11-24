function Y = NewMexico_calibrate_readoutput(output_dir,P)

switch P.FixedParameters.QoI

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
        r_exp = P.FixedParameters.r_exp; %2.25*[0.25 0.35 0.6 0.82 0.92]; % Measurement radial stations in percentage of blade length
        
        switch P.FixedParameters.QoI_type
            
            case 'mean'
                % Calculate mean
                Fn_mean = mean(Fn,1); % Mean (average in time) values at different radial stations

                % Interpolation
                Y   = spline(r_sim,Fn_mean,r_exp); % Interpolated data using spline                
                
            case 'full'
                % Use full (azimuth dependent) solution
                % select only the last revolution
                delta_azi = 10;
                azi_sim   = D{:,2};
                ind_last_rev = find(azi_sim<delta_azi,2,'last');
                Fn_last_rev = Fn(ind_last_rev(1):ind_last_rev(2)-1,:);
                azi_last_rev = azi_sim(ind_last_rev(1):ind_last_rev(2)-1);
                
                % Interpolation: columns of Fn_last_rev are interpolated to
                % yield new columns at r_exp positions
                Y   = spline(r_sim,Fn_last_rev,r_exp);

                % now interpolate to the azimuth positions of the
                % experimental data
                % use transpose to make interpolation of entire matrix
                % possible
                azi_exp = P.FixedParameters.azi_exp; %
                Y   = spline(azi_last_rev',Y',azi_exp)';
                
                % alternatively, we could directly do 2D interpolation
                
                % put all into a single row vector
                Y   = Y(:)';                
                
            otherwise
                error('QoI type unknown');
        end
        
    otherwise
        error(strcat('QoI type unknown'));
        
end
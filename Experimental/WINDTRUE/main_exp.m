close all
clear all
clc
format longG

% the .tim files contain time series: 133 sensors at 35 samples/sec
% with 10 minutes this gives 21000 points per run
dirname = ['../InflowLoadsAerodynamicsMeasurements_35Hz'];

listing = dir(dirname);

% loop over directory
for i = 1:size(listing,1)
    
    
    if (~strcmp(listing(i).name(1),'.')) % ignore files starting with .
        
        [filepath,namefile,extfile] = fileparts(listing(i).name);
        
        day = namefile;
        
        dirnameday = [dirname filesep day];
        
        % loop over subdirectories
        listingday = dir(dirnameday);
        
        for j = 1:size(listingday,1)
            
            [filepath,namefile,extfile] = fileparts(listingday(j).name);
            
            if (~strcmp(listingday(j).name(1),'.')) % ignore files starting with .
                
                ts = namefile;
                
                dirnamets = fullfile(dirname,day,[ts extfile]);
                disp(dirnamets)
                
                daqwin_noc = 174;  % Number of channels (see Table 9-1)
                
                % read the raw (binary) data:
                daqwin     = ReadData(dirnamets,daqwin_noc);
                
                % process some of the data
                daqwin(:,175) = daqwin(:,70)-daqwin(:,131); % Yaw error
                daqwin(:,138) = -daqwin(:,138); % Fx03
                daqwin(:,139) = -daqwin(:,139); % Fx05
                daqwin(:,140) = -daqwin(:,140); % Fx08
                daqwin(:,141) = -daqwin(:,141); % Fx10
                
                %         index = findconstrpm(daqwin,dirnamets);
                
                %         [daqwin] = filtconstrpm(daqwin,index);
                
                % select a certain time range, see also paper of
                % Madsen et al., Measured aerodynamic forces on a full scale 2MW turbine in comparison with
                % EllipSys3D and HAWC2 simulations
                [daqwin] = filt(daqwin,200,450);
                
                % write statistics
                statis(daqwin)
                
                %         axindstudy(daqwin)
                
                %         plaw(daqwin	
                
                [mFx,sFx] = plotFx(daqwin);
                
                [mFy,sFy] = plotFy(daqwin);
                
                vector = [mFy sFy mFx sFx];
                
                save('vector.mat','vector');
            end
            
        end
        
    end
end
%% Read the raw data

R1 = daqwin(:,142); % Fy03
R2 = daqwin(:,143); % Fy05
R3 = daqwin(:,144); % Fy08
R4 = daqwin(:,145); % Fy10
R5 = daqwin(:,138); % Fx03
R6 = daqwin(:,139); % Fx05
R7 = daqwin(:,140); % Fx08
R8 = daqwin(:,141); % Fx10

%%
N = table(mFy',sFy', 'VariableNames',{'mean','sd'});
writetable(N,'output_exp.dat','Delimiter','\t','WriteRowNames',true);
type output_exp.dat
T = table(mFy', 'VariableNames',{'exp_data'});
writetable(T,'output_e.dat','Delimiter','\t','WriteRowNames',true);
type output_e.dat
Raw = table(R1,R2,R3,R4, 'VariableNames',{'Fy03','Fy05','Fy08','Fy10'});
writetable(Raw,'raw.dat','Delimiter','\t','WriteRowNames',true);
type raw.dat;
Raw2 = table(R5,R6,R7,R8, 'VariableNames',{'Fx03','Fx05','Fx08','Fx10'});
writetable(Raw2,'raw_tan.dat','Delimiter','\t','WriteRowNames',true);
type raw_tan.dat;

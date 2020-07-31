close all
clear all
clc
format longG

dirname = ['../InflowLoadsAerodynamicsMeasurements_35Hz'];

listing = dir(dirname);

for i = 3:size(listing,1)
    
    [filepath,namefile,extfile] = fileparts(listing(i).name);
    
    day = namefile;
    
    dirnameday = [dirname filesep day];
    
    listingday = dir(dirnameday);
    
    for i = 3:size(listingday,1)
        
        [filepath,namefile,extfile] = fileparts(listingday(i).name);
        
        ts = namefile;
        
        dirnamets = [dirname filesep day filesep ts extfile]
        
        daqwin_noc=174;  % Number of channels (see Table 9-1)
        daqwin=ReadData(dirnamets,daqwin_noc);
        
        daqwin(:,175) = daqwin(:,70)-daqwin(:,131); % Yaw error
        daqwin(:,138) = -daqwin(:,138); % Fx03
        daqwin(:,139) = -daqwin(:,139); % Fx05
        daqwin(:,140) = -daqwin(:,140); % Fx08
        daqwin(:,141) = -daqwin(:,141); % Fx10
        
%         index = findconstrpm(daqwin,dirnamets);
        
%         [daqwin] = filtconstrpm(daqwin,index);
        
        [daqwin] = filt(daqwin);
        
        statis(daqwin)
        
%         axindstudy(daqwin)
        
%         plaw(daqwin)
        
        [mFx,sFx] = plotFx(daqwin);
        
        [mFy,sFy] = plotFy(daqwin);
        
        vector = [mFy sFy mFx sFx];
        
        save('vector.mat','vector');
        
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

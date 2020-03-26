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

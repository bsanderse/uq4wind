function [polar, Re] = chav_hansen_3D(k,x,section,c,r,pitch,twist)
    
    disp('CORRECTION STARTED!!!!')
    
    %% data input - read file
   
    
    fid = fopen(section,'r');
    
    if fid == -1 
        
        error('Error, check 2D polar file name')
        
    end
        
    for j = 1:8
                
        tline = fgetl(fid);
                
        if j == 8 
           
        Reyn = tline;
                   
        Re = str2double(extractAfter(Reyn,' '));
                   
        end
    end
       
    data = textscan(fid,'%f %f %f %f %f %f %f','Delimiter','tab');
    
    fclose(fid);
    
    aoa = [data{1}];
    
    cl2d_visc = [data{2}];
    
    cd2d = [data{3}];
    
    %% Find Cl_0 
     
    cl_0 = findzero(aoa, cl2d_visc);
    
    cl2d_inv = 2*pi*deg2rad(aoa) + cl_0;     
    
    %% Dcl, Dcd, f coefficient calculations 
    
    cd2d_min = minpositive(cd2d); %Find minimum 2D Cd (inviscid), minimum drag value
    
    Dcl = cl2d_inv - cl2d_visc;
    
    Dcd = cd2d - cd2d_min;
    
    theta = deg2rad(pitch + twist);
    
    f = k*(c/r)*(cos(theta))^x;
    
    %% create the 3D polar for each aoa
    
    rows = size(aoa,1);
    
    cl3d = zeros(rows,1);
    cd3d = zeros(rows,1);
      
    cl3d = cl2d_visc + f*Dcl;
    
    cd3d = cd2d + f*Dcd;
        
    polar = [aoa,cl3d,cd3d,data{4}];
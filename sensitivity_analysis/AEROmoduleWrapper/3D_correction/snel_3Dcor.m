function [polarm,Re] = snel_3Dcor(k,x,section,c,r)

    %%%% Snel model for 3D correction takes coefficients and creates 3D
    %%%% polar files in PhatAero_BEM
    
    % add the path to 2D polars
    addpath('3D_correction\2D_polars');
    
    %% data input - read file
    
    fid = fopen(section,'r');
    
    if fid == -1 
        
        disp('Error, check 2D polar file name')
        
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
    
    cd(oldFolder);
    
    aoa = [data{1}];
    
    cl2d_visc = [data{2}];
    
    %% Find Cl_0
     
    cl_0 = findzero(aoa, cl2d_visc);

    
    cl2d_inv = 2*pi*deg2rad(aoa) + cl_0; 
    
    
    %% Viscous - Inviscid difference
    
    Dcl = cl2d_inv - cl2d_visc;
    
    
        %% 3D Cl calculation
        
        rows = size(aoa,1);
        fcl = zeros(rows,1); %add it with size
        cl3d = zeros(rows,1);
        
        
            
        %%% Calculation of the coefficient
    
        fcl = k*(c/r)^x; 
    
        %%% 3D correction
    
        cl3d = cl2d_visc + fcl*Dcl;
        
        polarm = [aoa,cl3d,data{3},data{4}];
        
    end
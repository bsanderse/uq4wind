 function write_polars(type,kfactor,exp,pitch)

%%% Chaviaropoulos - Hansen model currently only works for constant pitch

%% Temporary bug fix with paths wr_path points always inside sensitivity_analysis folder

 root_folder = pwd;
 
 if endsWith(root_folder,'NewMexico_calibrate_allruns') 
     
    wr_path = fullfile(root_folder, '..\..');
   
    
 elseif endsWith(root_folder, 'sensitivity_analysis')
     
    wr_path = root_folder;
    
 else 
     error('Error check paths that function write_polars uses')
     
 end
     
     
%% 3D Snel, type == 1
    
   if type == 1
        
        dinfo = dir(fullfile(wr_path,'AEROmoduleWrapper\3D_correction\2D_polars'));
        
        files = {dinfo.name};
        
        files1 = {files{1,3:7}}; %extract the 2D polars filenames
        
        fid = fopen('airfoils_g.dat','r');

        if fid == -1 

            disp('Error, check file name')

        end
        
        headers = fgetl(fid);
        
        geom = textscan (fid,'%s %f %f %f %f','Delimiter','\t');
        
        section = geom{1,1};
        r = geom{1,2};
        tpc = geom{1,3};
        c = geom{1,4};
        
        %%% Model coefficients
        k = kfactor;
        x = exp;

        for i = 1 : size(files1,2)

            [polar_mat, Re(i)] = snel_3Dcor(k,x,files1{i},c(i),r(i));
           
            %%%%%check if folder changes
            
            filevar = fullfile(wr_path,'AEROmodule\NewMexico_calibrate\current',section{i}); % path to polar in current folder for Aeromodule
            
            fileID = fopen(filevar,'w');

            if fileID == -1

                disp('Error, check 3D polar file name')

            end

            %%%%% Write the required format for polars
            fprintf(fileID,'! Aero module input file for airfoil data\n');
            fprintf(fileID,'\n');
            fprintf(fileID,'Airfoil_Name  %s\n',section{i});
            fprintf(fileID,'t/c           %s         ! thickness ratio w.r.t. chord\n',tpc(i));
            fprintf(fileID,'\n');
            fprintf(fileID,'format 1       !  1: alfa-cl-cd-cm	; 2: alfa-cl; alfa-cd; alfa-cm\n');
            fprintf(fileID,'\n');
            fprintf(fileID,'Reynolds_Nr %f\n',Re(i));
    
            writematrix(polar_mat, filevar, 'WriteMode', 'Append', 'Delimiter', '\t');
            
            disp([section{i} ' polar file created!'])
            
            fclose(fileID);
            

        end  
    
    %% Chaviaropoulos - Hansen 3D, type == 2
    
    elseif type == 2
        
        dinfo = dir(fullfile(wr_path,'AEROmoduleWrapper\3D_correction\2D_polars'));
        
        files = {dinfo.name};
        
        files1 = {files{1,3:7}}; %extract the 2D polars
        
        
        fid = fopen('airfoils_g.dat','r');
        
        if fid == -1

            disp('Error, check file name')

        end
        
        headers = fgetl(fid);
        
        geom = textscan (fid,'%s %f %f %f %f','Delimiter','\t');
        
        section = geom{1,1};
        r = geom{1,2};
        tpc = geom{1,3};
        c = geom{1,4};
        twist =geom{1,5};

        % 3D correction model coefficients

        k = kfactor;
        x = exp;

        %for each section
        for i = 1 : size(files1,2)

            [polar_mat, Re(i)] = chav_hansen_3D(k,x,files1{i},c(i),r(i),pitch,twist(i));

            %%%%%check if folder changes
            filevar = fullfile(wr_path,'AEROmodule\NewMexico_calibrate\current',section{i}); % path to polar in current folder for Aeromodule
            
            fileID = fopen(filevar,'w');

            if fileID == -1

                disp('Error, check 3D polar file name')

            end

            %%%%% Write the required format for polars
            fprintf(fileID,'! Aero module input file for airfoil data\n');
            fprintf(fileID,'\n');
            fprintf(fileID,'Airfoil_Name  %s\n',section{i});
            fprintf(fileID,'t/c           %s         ! thickness ratio w.r.t. chord\n',tpc(i));
            fprintf(fileID,'\n');
            fprintf(fileID,'format 1       !  1: alfa-cl-cd-cm	; 2: alfa-cl; alfa-cd; alfa-cm\n');
            fprintf(fileID,'\n');
            fprintf(fileID,'Reynolds_Nr %f\n',Re(i));
            
            writematrix(polar_mat, filevar, 'WriteMode', 'Append', 'Delimiter', '\t');
            
            disp([section{i} ' polar file created!'])
            
            fclose(fileID);

        end

        
    
    
    %% Error, wrong type
    
    else
        
        error('Error, invalid model type')
    
    end
    
end
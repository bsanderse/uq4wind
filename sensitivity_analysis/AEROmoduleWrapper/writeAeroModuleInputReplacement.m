function writeAeroModuleInputReplacement(X,P)
% This routine adds the random input to the input.txt file used for Aero
% module.

FixedParameters = P.FixedParameters;
UncertainInputs = P.UncertainInputs;

root_folder    = FixedParameters.root_folder;
ref_folder     = FixedParameters.ref_folder;
current_folder = FixedParameters.current_folder;


%% copy all files from reference folder to current folder
ref_dir = fullfile(root_folder,ref_folder);
cur_dir = fullfile(root_folder,current_folder);
copyfile(fullfile(ref_dir,'*.dat'),cur_dir);
copyfile(fullfile(ref_dir,'*.txt'),cur_dir);
copyfile(fullfile(ref_dir,'*.ini'),cur_dir);


%% now adapt files that need to be changed due to uncertainties
%loop over uncertain inputs and replace the ones in input.txt

% open the reference input.txt file
filename_in  = fullfile(cur_dir,'input.txt');
fid_in       = fopen(filename_in,'r');
lines        = textscan(fid_in,'%s','delimiter','\n');
fid_in       = fclose(fid_in);

lines = lines{1};
lines_new = lines;

zz = 1;

write_wind_file = 0; % determines whether the wind.dat file needs to be updated

% loop over all uncertainties (and constants), and look for a match in the
% input files
% keep track if all uncertainties are used
ndim = length(X);
uncertainty_covered = zeros(ndim,1);

for i=1:ndim
    UncertainInputName = UncertainInputs.Marginals(i).Name;
    
    switch UncertainInputName
        
        case {'CL','CD','CM'}
            
            %% update the airfoil polar files
            % administrate which polars are written
            index_pol       = UncertainInputs.Marginals(i).AirfoilIndex;
            polar_index(zz) = index_pol;
            
            % airfoil filename
            airfoil{index_pol}  = UncertainInputs.Marginals(i).Airfoil;
            airfoil_file = fullfile(cur_dir,airfoil{index_pol});
            
            % perturbed indices
            alpha_pert = UncertainInputs.Marginals(i).AlphaPert;
            
            % check if this section was already loaded before
            % if so, then don't load again to prevent overwriting the
            % perturbed polars
            if (length(polar_index)==length(unique(polar_index)))
                [alpha{index_pol}, CL{index_pol}, CD{index_pol}, CM{index_pol}, header{index_pol}] = importPolars(airfoil_file);
            end
            
            % computeCurves constructs a NURBS in Cl-alpha space
            % alpha_pert consists of the indices that are perturbed
            ind_alpha_pert = find(alpha{index_pol}>=alpha_pert(1) & alpha{index_pol}<=alpha_pert(2));
            d = ones(length(ind_alpha_pert),1);
            plotCurve = 0;
                        
            switch UncertainInputName
                
                case 'CL'
                    CL_pert = computeCurves(1, ind_alpha_pert, X(index_pol)*d, d, plotCurve, ...
                        alpha{index_pol}, CL{index_pol}, 3, 1:length(alpha{index_pol}));                    
%                     figure(i+100)
%                     hold on
%                     plot(alpha{index_pol},CL_pert,'--','LineWidth',2);
                    
                    CL{index_pol}   = CL_pert;
                    
                case 'CM'
                    
                    CM_pert = computeCurves(1, ind_alpha_pert, X(index_pol)*d, d, plotCurve, ...
                        alpha{index_pol}, CM{index_pol}, 3, 1:length(alpha{index_pol}));
                    
%                     figure(i+100)
%                     hold on
%                     plot(alpha{index_pol},CM_pert,'--','LineWidth',2);                    
                    
                    CM{index_pol}   = CM_pert;

                case 'CD'
                    
                    CD_pert = computeCurves(1, ind_alpha_pert, X(index_pol)*d, d, plotCurve, ...
                        alpha{index_pol}, CD{index_pol}, 3, 1:length(alpha{index_pol}));
                    
%                     figure(i+100)
%                     hold on
%                     plot(alpha{index_pol},CD_pert,'--','LineWidth',2);                    

                    CD{index_pol}   = CD_pert;

            end
            zz = zz + 1;
            uncertainty_covered(i) = 1;
            
        case {'Twist','Chord','Thickness'}
            
            %% update the AEROPROPS in the input.txt file
            ind1 = find(startsWith(lines,'AEROPROPS'));
            if (~isempty(ind1) && ind1>0)
                % try to determine how many lines the AEROPROPS table has
                i1 = find(startsWith(lines,'BLADEROOT'));
                i2 = find(startsWith(lines,'BLADELENGTH'));
                i3 = find(startsWith(lines,'AEROROOT'));
                ind2 = min([i1 i2 i3]);
                if (strcmp(lines{ind1+1}(1),'!')) % skip over the !zB line
                    ind1 = ind1+1;
                end
                if (find(startsWith(lines(ind1+1:ind2-1),'!')))
                    error('please cleanup the AEROPROPS table');
                end
                % we now know that we should replace lines ind1+1:ind2-1
                ncol  = length(str2num(lines{ind1+1})); % generally we have 7 columns
                nrow  = (ind2-1)-(ind1+1)+1;
                aeroprops = zeros( nrow, ncol);
                for k=1:nrow
                    aeroprops(k,:) = str2num(lines{k + ind1});
                end
                
                % select Twist, Chord, or Thickness
                
                error('code to be finished');
                
                % sprintf(fid,'%f    %f    %f    %f    %f    %f    %f \n', P{3}(i), chord(i), thickness(i)/chord(i), twist(i), P{7}(i), P{8}(i), P{9}(i));
                
            else
                error('no AEROPROPS keyword in input.txt');
            end
            uncertainty_covered(i) = 1;
            
        case 'WINDSPEED'
            write_wind_file = 1;
            wind_file_line_id = find(startsWith(lines,'WINDFILENAME'));
            wind_file_line = strsplit(lines{wind_file_line_id}); % split line in 2; the WINDFILENAME and the actual  name
            wind_file = wind_file_line{2}; % second entry on the line
            V_inf = X(i);
            uncertainty_covered(i) = 1;
            
        case {'factor3D', 'exp3D'}
            
           %% case for testing Chaviaropoulos - Hansen model 
            factor3D = X(1);
            exp3D = X(2);
            
            % Choose 3D correction type, 1 --> Snel, 2 --> Chav.-Hansen
            
            type = 2;
            
            if type == 1
                
                disp('Snel correction selected')
                
            elseif type == 2
                
                disp('Chaviaropoulos - Hansen correction selected')
            
            else
                
                disp('Invalid type of correction')
            end
            
            write_polars(type,factor3D,exp3D)
            
            uncertainty_covered(i) = 1;
            

            
        otherwise
            %% generic scalar variables
            % check if a line starts with the variable name (case sensitive)
            % this will skip any commented lines, i.e. those that start with !
            
            ind = find(startsWith(lines,UncertainInputName));
            if (~isempty(ind) && ind>0)
                lines_new{ind} = [UncertainInputName '      ' num2str(X(i))];
                uncertainty_covered(i) = 1;
            end
            
    end
    
    
end

% write new lines to input.txt in current folder
filename_out = fullfile(cur_dir,'input.txt');
fid_out      = fopen(filename_out,'w');
for i = 1:length(lines)
    fprintf(fid_out,'%s\n',lines_new{i});
end
fid_out = fclose(fid_out);

% now write the polar files
if (exist('polar_index','var'))
    polar_index = unique(polar_index);
    for q = 1:length(polar_index)
        i = polar_index(q);
        airfoil_out = fullfile(cur_dir,airfoil{i});
        % open the polar file
        fid_polar = fopen(airfoil_out,'w');
        % wrtie header
        for j = 1:length(header{i})
            fprintf(fid_polar,'%s \n',header{i}{j});
        end
        % print polar
        for j = 1:length(alpha{i})
            fprintf(fid_polar,'%f    %f    %f    %f \n', alpha{i}(j), CL{i}(j), CD{i}(j), CM{i}(j));
        end
        fclose(fid_polar);
    end
end


%% repeat but now for specialist_input.txt

% check if there is a specialist_input file
specialist_file_line_id = find(startsWith(lines,'INCLUDE'));
if (specialist_file_line_id>0)
    specialist_file_line = strsplit(lines{specialist_file_line_id}); % split line in 2; the WINDFILENAME and the actual  name
    specialist_file = specialist_file_line{2}; % second entry on the line
    
    % open the reference specialist_input.txt file
    filename_in  = fullfile(cur_dir,specialist_file);
    fid_in       = fopen(filename_in,'r');
    
    if (fid_in>0) % check if file exists / was opened successfully
        lines        = textscan(fid_in,'%s','delimiter','\n');
        fid_in       = fclose(fid_in);
        
        lines = lines{1};
        lines_new = lines;
        
        for i=1:ndim
            % check if a line starts with the variable name (case sensitive)
            % this will skip any commented lines, i.e. those that start with !
            UncertainInputName = UncertainInputs.Marginals(i).Name;
            
            % in this case we don't look at polars or AEROPROPS but only at
            % generic scalar variables
            ind = find(startsWith(lines,UncertainInputName));
            if (~isempty(ind) && ind>0)
                lines_new{ind} = [UncertainInputName '      ' num2str(X(i))];
                uncertainty_covered(i) = 1;
            end
            
        end
        
        % write new lines to specialist_input.txt
        filename_out = fullfile(cur_dir,'specialist_input.txt');
        fid_out      = fopen(filename_out,'w');
        for i = 1:length(lines)
            fprintf(fid_out,'%s\n',lines_new{i});
        end
        fid_out = fclose(fid_out);
    else
        warning('note: specialist input file could not be opened');
    end
else
    disp('note: no specialist input file');
end


%% write the wind.dat file if necessary
if (write_wind_file == 1)
    
    % open the wind.dat file
    filename_wind  = fullfile(cur_dir,wind_file);
    
    % read table
    wind_data = readtable(filename_wind);
    % first column: time, second column: u, third: v, fourth: w
    % set second column equal to specified value
    wind_data(:,2) = table(V_inf*ones(height(wind_data),1));
    
    % write to table, note that headers are skipped
    writetable(wind_data,filename_wind,'Delimiter','space','WriteVariableNames',false);
    
end

%% check if there are any uncertainties not written to one of the files
ind_notcovered = find(uncertainty_covered == 0);
if (~isempty(ind_notcovered))    
    for i=1:length(ind_notcovered)
        warning([UncertainInputs.Marginals(ind_notcovered(i)).Name ' not used when writing input files']);
    end    
end

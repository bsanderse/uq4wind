%% Script to extract data from NewMexico force measurements in the required format of the WINDTRUE framework
   % the files have to be in the format originally supplied as of(29/01/21) by TNO through
   % the folder \\tsn.tno.nl\RA-Data\SV\sv-092280\DEMO\ECN_reduced_data\Bin Average Forces
   % Run this script inside the folder with the required data
%% Clear
clc
clearvars
close all

root_folder = pwd;
%% Get filenames of csv saved force measurements
d = dir('*.csv');
filenames = {};
index = 0;
k = 0;
for i = 1:length(d)
  % Using regexp with lookaround on the file name to find 
  % a sequence of digits preceeded by 'D' and followed by '_loads': 
  no = regexp(d(i).name,'(?<=D)\d*(?=\_loads)','match');
  if ~isempty(no)
    k = k+1;
    filenames{k} = d(i).name;
    index(k) = str2num(no{end});
  end
end

%% Create new txt files with the desired format for the framework

num_files = length(filenames);

for i = 1:num_files
    
    %% Open file
    fid = fopen(filenames{i},'r');
    
    if fid == -1 

            disp('Error, check file name')

    end
    %% Extract header lines and remove std and avg lines

                
    tline = fgetl(fid);
       
    headers = tline;
    
    data = readtable(filenames{i}, 'HeaderLines', 3); %skip first 3 rows
    
    fclose(fid);
    
    %% Create a txt file with the desired format for each csv 
    
    filevar = strcat(erase(filenames{i},'.csv'), '.dat');        
    
    fid = fopen(filevar,'w');
    
    % remove unused column headers
    headers = erase(headers, 'Fa');
    headers = erase(headers, 'Torque');
    
    fprintf(fid, headers);
    
    datamat = data(:,[1 4:13]);
  
    writetable(datamat, filevar, 'WriteMode', 'Append', 'Delimiter', '\t');
    
    fclose(fid);
end   
    
    
    
    
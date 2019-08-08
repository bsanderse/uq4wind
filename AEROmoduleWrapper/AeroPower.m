function [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [TIMES,AZIMUTHDEG,POWERWATT,AXIAL_FORCEN] = AeroPower(FILENAME) Reads
%   data from text file FILENAME for the default selection.
%
%   [TIMES,AZIMUTHDEG,POWERWATT,AXIAL_FORCEN] = AeroPower(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower('AeroPower.dat',2, 445);

%% Initialize variables.
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%20f%20f%20f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Times = dataArray{:, 1};
Azimuthdeg = dataArray{:, 2};
PowerWatt = dataArray{:, 3};
Axial_ForceN = dataArray{:, 4};



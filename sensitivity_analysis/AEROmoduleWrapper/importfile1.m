function outpute = importfile1(filename, dataLines)
%IMPORTFILE1 Import data from a text file
%  OUTPUTE = IMPORTFILE1(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  OUTPUTE = IMPORTFILE1(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  outpute = importfile1("C:\Users\kriek\Desktop\windtrue\Experimental\WINDTRUE\output_e.dat", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 22-Apr-2020 15:28:30

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "exp_data";
opts.VariableTypes = "double";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
outpute = readtable(filename, opts);

end
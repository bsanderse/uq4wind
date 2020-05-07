function outputexp = importfile2(filename, dataLines)
%IMPORTFILE2 Import data from a text file
%  OUTPUTEXP = IMPORTFILE2(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  OUTPUTEXP = IMPORTFILE2(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  outputexp = importfile2("C:\Users\kriek\Desktop\windtrue\Experimental\WINDTRUE\output_exp.dat", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 06-May-2020 14:44:39

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["mean", "sd"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
outputexp = readtable(filename, opts);

end
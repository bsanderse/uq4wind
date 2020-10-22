function loads = readNewMexico(filename, dataLines)
%READNEWMEXICO Import data from a text file
%  LOADS = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  LOADS = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  loads = importfile("\\mac\home\Dropbox\work\Programming\UQ\windtrue\Experimental\NewMexicoData\R52P78D929_loads.dat", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 21-Oct-2020 09:21:21

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Azi", "Fn25Npm", "Fn35Npm", "Fn60Npm", "Fn82Npm", "Fn92Npm", "Ft25Npm", "Ft35Npm", "Ft60Npm", "Ft82Npm", "Ft92Npm"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
loads = readtable(filename, opts);

end
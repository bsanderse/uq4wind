function [alpha, Cl, Cd, Cm] = importPolars(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [ALPHA, CL, CD, CM] = importPolars(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as column
%  vectors.
%
%  [ALPHA, CL, CD, CM] = importPolars(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  [alpha, Cl, Cd, Cm] = importPolars("..\AEROmodule\NM80\section03.dat", [9, Inf]);

% %% Input handling
% 
% % If dataLines is not specified, define defaults
% if nargin < 2
%     dataLines = [9, Inf];
% end
% 
% %% Setup the Import Options
% 
% if (verLessThan('matlab','9.5'))
%     % for older versions of Matlab:
%     tbl = readtable(filename, 'Delimiter', '\t');
%     tbl.Properties.VariableNames = {'alpha','Cl','Cd','Cm'};
%     
% else
%     
%     % opts = delimitedTextImportOptions("NumVariables", 4);
%     
%     % Specify range and delimiter
%     opts.DataLines = dataLines;
%     opts.Delimiter = "\t";
%     
%     % Specify column names and types
%     opts.VariableNames = ["alpha", "Cl", "Cd", "Cm"];
%     opts.VariableTypes = ["double", "double", "double", "double"];
%     opts.ExtraColumnsRule = "ignore";
%     opts.EmptyLineRule = "read";
%     
%     % Import the data
%     pwd
%     disp(filename)
%     tbl = readtable(filename, opts);
%     
% end
% 
% %% Convert to output type
% alpha = tbl.alpha';
% Cl = tbl.Cl';
% Cd = tbl.Cd';
% Cm = tbl.Cm';
% 
% %% Import data from text file
% % Script for importing data from the following text file:
% %
% %    filename: C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80\section03_ref.dat
% %
% Auto-generated by MATLAB on 15-Apr-2020 16:50:39

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [9, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["alpha", "Cl", "Cd", "Cm"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
alpha = tbl.alpha';
Cl = tbl.Cl';
Cd = tbl.Cd';
Cm = tbl.Cm';

 


%% Clear temporary variables
clear opts
end
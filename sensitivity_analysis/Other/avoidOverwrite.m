function outFile = avoidOverwrite(inFile,inPath,numDigits,startCount)
% Function to check if a given file already exists and create a new
% filename that can be used to avoid overwriting. New file names are
% created by appending a sequence like '_001' to the file name. The number
% is hereby increased until an unused file name is found. 
%
% Usage: outFile = avoidOverwrite(inFile,inPath,numDigits)
%
%        inFile: The name of the file that should be checked. Make sure
%        to add filetype if required.
%        inPath: The path of the file. Can be omitted to only check in
%        the Matlab path.
%        numDigits: The number of digits used for enumration. Default is 2
%        digits, so filenames are set up as e.g. 'inFile_01.mat' and so on.
%        startCount: Defines the starting number for enumeration. Standard
%        is 0 so the first file is created as 'inFile_00.mat'. 
%        outFile: The filename that can be used to save a file without
%        overwriting existing data.
%       
%        Example: 
%        fName = 'Example.mat'; %example file
%        outFile1 = avoidOverwrite(fName,pwd,3,1);
%        save(outFile1); %save empty file
%        disp(['First file name: ' outFile1]);
%        outFile2 = avoidOverwrite(fName,pwd,3,1); %second file
%        save(outFile2);
%        disp(['Second file name: ' outFile2]);
%        outFile3 = avoidOverwrite(fName,pwd,3,1); %third file
%        save(outFile3);
%        disp(['Third file name: ' outFile3]);
%        delete(outFile1,outFile2,outFile3) %delete empty files again

%% check input
if ~exist('inPath','var')
    inPath = [];     %number of digits when adding numbers to a filename
end

if ~exist('numDigits','var')
    numDigits = 2;     %number of digits when adding numbers to a filename
end
numDigits = num2str(numDigits);

if ~exist('startCount','var')
    startCount = 0;     %first value for counter. This determines the first filename during enumeration.
end

if ~strcmp(inPath(end),filesep) %make sure path has a seperator at the end
    inPath = [inPath filesep];
end

%% check if file exists already and enumerate
if exist([inPath inFile],'file') == 2 %file exists already, check for alternative
    [~,shortFile,fileType] = fileparts(inFile); %exclude the file type
    checker = true; %check for alternate file names
    Cnt = startCount; %counter for file name
    
    while checker
        testPath = [inPath shortFile '_' num2str(Cnt, ['%0' numDigits 'i']) fileType];
        
        if exist(testPath,'file') == 2
            Cnt = Cnt + 1; %increase counter until a non-existing file name is found
        else
            checker = false;
        end
        
        if Cnt == 10^numDigits-1 && checker
            numDigits = numDigits+1;
            warning(['No unused file found at given number of digits. Number of digits increased to ' num2str(numDigits) '.']);
        end
    end
    outFile = [shortFile '_' num2str(Cnt, ['%0' numDigits 'i']) fileType];
    
else
    outFile = inFile;
end
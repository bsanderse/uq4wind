function path_found = findAeroModulePath()
path_strings = strsplit(getenv('PATH'),';');
path_found   = false;
for i=1:length(path_strings)
    if (strfind(path_strings{i},'ECNAero')>0)
        disp(['ECNAero location: ' path_strings{i}]);
        path_found = true;
    end
end
if (~path_found)
    error('ECNAero path not found');
end
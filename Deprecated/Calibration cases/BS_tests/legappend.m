function legappend(newStrings)
% Quick & dirty fork of legappend specific to MATLAB versions >= R2014b
% Only supports appending strings to the existing legend handle
% Assumes only one legend is in the current figure
% Add multiple strings by passing it a 1D cell array of strings

% Find our legend object
h = findobj(gcf, 'Type', 'legend');

if ~isempty(h)
    % Get existing array of legend strings and append our new strings to it
    oldstr = h.String;
    if ischar(newStrings)
        % Input string is a character array, assume it's a single string and
        % dump into a cell
        newStrings = {newStrings};
    end

    newstr = [oldstr{1:end-1} newStrings]; % BS 2020: remove the automatic entry from Matlab

    % Get line object handles
    ploth = flipud(get(gca, 'Children'));

    % Update legend with line object handles & new string array
    h.PlotChildren = ploth;
    h.String = newstr;
end
end
% 
% function [legend_h,object_h,plot_h,text_strings] = legappend(newStrings,varargin)
% %LEGAPPEND appends new entries to the end of a legend by deleting the
% %current legend and recreating a new, similar legend. 
% % 
% %% Syntax
% % 
% % legappend('new legend entry') 
% % legappend('new entry 1','new entry 2',...,'new entry N') 
% % legappend('') 
% % legappend('','',...,'')
% % [legend_h,object_h,plot_h,text_strings] = legappend(...)
% % 
% %% Description 
% % 
% % legappend('new legend entry') appends an existing legend with "new
% % legend entry".
% % 
% % legappend('new entry 1','new entry 2',...,'new entry N') adds several
% % new entries to the legend. 
% % 
% % legappend('') deletes the last entry from the legend. 
% % 
% % legappend('','',...,'') deletes the last several entries from the
% % legend.
% % 
% % [legend_h,object_h,plot_h,text_strings] = legappend(...) returns legend_h, the
% % handle of the new legend; object_h, handles of the line, patch, and
% % text graphics objects used in the legend; plot_h, handles of the lines
% % and other objects used in the plot; and text_strings, a cell array of
% % the text strings used in the legend. Note that for new legend entries,
% % legappend does not add entries to a current legend, but deletes the
% % current legend and recreates a new one. As a result, the legend handle
% % will change with each new-entry use of legappend.  The legend handle
% % does not change when legappend is used to delete an entry. 
% % 
% % 
% %% Author Info
% % This function was created by Chad A. Greene of the Institute for
% % Geophysics, The University of Texas at Austin, July 2014. 
% % 
% % See also legend.
% 
% % BS 2018: change for newer Matlab version as suggested on File Exchange
% % h =  findobj(gcf,'Type','axes','Tag','legend');
% 
% h = findobj(gcf, 'Type', 'legend'); 
% newstr = [h.String {newStrings}]; 
% allDatah = flipud(get(gca,'children')); 
% h.PlotChildren = allDatah; 
% h.String = newstr;
% 
% prop.boxon = get(h,'visible');
% prop.loc = get(h,'location'); 
% prop.color = get(h,'color'); 
% prop.orient = get(h,'Orientation'); 
% 
% 
% % allDatah = flipud(get(gca,'children')); 
% str = get(h,'String'); 
% 
% if exist('varargin','var') 
%     newStrings = [newStrings,varargin];
% end
% deleteEntries = sum(cellfun('isempty',newStrings));
% if isempty(newStrings) 
%     deleteEntries = 1; 
% end
% 
% if ~deleteEntries
%     if iscell(newStrings)
%         for k = 1:length(newStrings) 
%             str{end+1}=newStrings{k}; 
%         end
%     end 
%     if ~iscell(newStrings)
%         str{end+1}=newStrings; 
%     end
% 
% 
%     [legend_h,object_h,plot_h,text_strings] = legend(h,allDatah,str);
% 
%     if strcmpi({prop.boxon},'off')
%         legend boxoff
%     end
% 
%     set(legend_h,'location',prop.loc,'color',prop.color,'Orientation',prop.orient)
% 
% 
% end
% 
% if deleteEntries
%     set(h,'String',str(1:end-nargin))
%     [legend_h,object_h,plot_h,text_strings] = legend;
% end
% 
% 
% if nargout==0
%     clear legend_h object_h plot_h text_strings
% end
% 
% 
% end


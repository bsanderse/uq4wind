function DefaultFont = uq_getDefaultFont(ax)
%UQ_GETDEFAULTFONT returns the default font properties of UQLab graphics.
%
%   DefaultFont = uq_getDefaultFont(AX) returns the default font properties
%   of UQLab graphics in a structure. For proper flexible sizing, the axes
%   object AX is required.
%
%   The following are the defaults:
%       FontSize                relative to figure size (0.04 * height,
%                               operation in centimeters)

%% Font properties

% FontSize property
% store settings we want to modify in here
fig = get(ax,'Parent');
currFigUnits = get(fig,'Units');
currAxFontUnits = get(ax,'FontUnits');
% Set the units to centimeters, scaling & adjustment are done in this unit
set(fig, 'Units', 'centimeters')
set(ax, 'FontUnits', 'centimeters')
sizeFactor = 0.04;  % ratio fontsize/figure-height
figSize = get(fig,'Position');  % size of current figure
% Set FontSize relative to figure size
fontSizeCm = sizeFactor*figSize(4);  % sizeFactor * Fig. height
fontSizeCurrUnits = convert_cmToX(fontSizeCm,currAxFontUnits);

DefaultFont.FontSize = fontSizeCurrUnits;

% Reset the modified 'Units' property
set(fig, 'Units', currFigUnits);
set(ax, 'FontUnits', currAxFontUnits);

end

%% ------------------------------------------------------------------------
function y = convert_cmToX(x,units)

inch2cm = 2.54;
inch2pt = 72;
inch2px = get_inch2px();
inch2nm = get_x2y('inches','normalized');

switch units
    case 'centimeters'
        y = x;
    case 'inches'
        y = x/inch2cm; 
    case 'points'
        y = x/inch2cm * inch2pt;
    case 'pixels'
        y = x/inch2cm * inch2px;
    case 'normalized'
        y = x/inch2cm * inch2nm;
end

end

%% ------------------------------------------------------------------------
function inch2px = get_inch2px()
% Get the conversion ratio from one inch to one pixel.

if ispc
    inch2px = 96;
end
if ismac
    inch2px = 72;
end
if isunix
    inch2px = get_x2y('inches','pixels');
end

end

%% ------------------------------------------------------------------------
function x2y = get_x2y(x,y)
%Get the conversion factor from 'x' unit to 'y' unit.

% Store current units
currUnits = get(0, 'Units');
% Set the units to pixels
set(0, 'Units', x)  
% Obtain this pixel information
scrSizeX = get(0,'ScreenSize');
% Sets the units of your root object (screen) to inches
set(0, 'Units', y)
% Obtains this inch information
scrSizeY = get(0,'ScreenSize');
% Calculates the resolution (pixels per inches)
x2y = scrSizeY(3)/scrSizeX(3);
% Restore the current units
set(0, 'Units', currUnits)

end

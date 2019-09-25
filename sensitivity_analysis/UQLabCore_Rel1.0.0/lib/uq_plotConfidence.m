function handles = uq_plotConfidence(X, Y, Conf, varargin)
% UQ_PLOTCONFIDENCE plots the confidence bounds: X versus Y with confidence
% CONF
% 
% See also: UQ_KRIGING_DISPLAY

% Set to default the coloring options:
Lcolor = uq_cmap(1);

%gray
BoundColor = [128,128,128]/255; 

Tag = '';
AlphaN = 0.3;
Lwidth = 2;
CreateFig = false;
if nargin > 3
    if mod(length(varargin), 2)
        error('The number of options is incorrect');
    end
    for ii = 1:2:length(varargin) - 1
        
        switch(lower(varargin{ii}))
            case 'linecolor'
                Lcolor = varargin{ii + 1};
            case {'tag', 'name'}
                Tag = varargin{ii + 1};
            case 'boundcolor'
                BoundColor = varargin{ii + 1};
            case 'alpha'
                AlphaN = varargin{ii + 1};
            case 'linewidth'
                Lwidth = varargin{ii + 1};
            case 'createfigure'
                CreateFig = varargin{ii + 1};
        end
    end
end
if isempty(Tag)
    MTag = '';
else
    MTag = ' - ';
end
if CreateFig
    uq_figure('Name', sprintf('%s%sConfidence plot', MTag, Tag));
end


hold on

% Plot the confidence bound:
handles(2) = fill([X(:); flipud(X(:))],[Y(:) + Conf(:); flipud(Y(:) - Conf(:))],'r','LineStyle','none');
handles(1) = plot(X, Y, 'Color', Lcolor, 'LineWidth', Lwidth);
set(handles(2), 'FaceColor', BoundColor)
alpha(AlphaN);
box on


% Remove the hold on:
hold off
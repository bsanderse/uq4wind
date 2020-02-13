function handle = uq_subplot(NRow,NCol,index)
%% UQ_SUBPLOT returns axes handle in tiled position
% specify margins
xMarginLab = 0.15; % label margin 
yMarginLab = 0.15; 
xMarginTot = 0.05; % total margin
yMarginTot = 0.05;
xMarginLoc = xMarginTot/(NCol);
yMarginLoc = yMarginTot/(NCol);
xLength = (1-xMarginLab*2-xMarginTot)/NCol;
yLength = (1-yMarginLab*2-yMarginTot)/NRow;

% x/y index
[xIndex,yIndex] = ind2sub([NRow,NCol],index);

% coordinates
xCoord = xMarginLab + xMarginLoc*(xIndex-1) + xLength*(xIndex-1);
yCoord = 1 - (yMarginLab + yMarginLoc*yIndex + yLength*(yIndex));

% create axes
handle = axes('Position',[xCoord yCoord xLength yLength]);
end

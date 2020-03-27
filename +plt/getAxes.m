function [axesHandle, figureHandle] = getAxes(axesOpt, idx)
%% A method for to create and/or select figure and axes for an upcoming plot based on some
% INPUTS(default values)
% axesOpt('res')--------String indicating the type of data to plot. Options are
%	.totalNumOfAxes(1)-----------------The total number of axes in the figure
%	.figureSize([400 400])-------------Size of the figure in pixels [width height]
%	.gapBetweenAxes(25)----------------Space (in pixels) between the axes
%	.btlrMargins([50 50 50 50])--------The total number of axes in the figure
%	.reposition(1)---------------------Whether to position figure at the top-right of the left-most screen
%	.numOfRows([])---------------------Specify the number of rows to use (for figures with multiple axes)

% OUTPUTS
% axesHandle-------------handle for genereated axis
% figureHandle-----------handle for the figure created/used

%% Set defaults etc.
%Set detault values
if ~exist('idx', 'var'); idx = 1; end
if ~isfield(axesOpt, 'totalNumOfAxes'); axesOpt.totalNumOfAxes = 1; end
if ~isfield(axesOpt, 'figureSize'); axesOpt.figureSize = [400 400]; end
if ~isfield(axesOpt, 'gapBetweenAxes'); axesOpt.gapBetweenAxes = 25; end
if ~isfield(axesOpt, 'btlrMargins'); axesOpt.btlrMargins = 50*ones(1,4); end
if ~isfield(axesOpt, 'reposition'); axesOpt.reposition = 1; end
if ~isfield(axesOpt, 'numOfRows'); axesOpt.numOfRows = []; end

%Detect the screensize and screen ratio (width/height)
screenSize = get(0,'MonitorPositions');
screenSize = screenSize(screenSize(:,1)==min(screenSize(:,1)),:);
screenRatio = round(screenSize(3)/screenSize(4));

if ~isempty(axesOpt.numOfRows); numOfRows = axesOpt.numOfRows;
elseif axesOpt.totalNumOfAxes < 4; numOfRows = 1;
else, numOfRows = find(4*screenRatio.*(1:100)>=axesOpt.totalNumOfAxes,1);
end

numOfCols = ceil(axesOpt.totalNumOfAxes/numOfRows);
axesOpt.figureSize = min([axesOpt.figureSize(1)*numOfCols, axesOpt.figureSize(2)*numOfRows], screenSize(3:4));

axesOpt.btlrMargins(1:2) = axesOpt.btlrMargins(1:2)/axesOpt.figureSize(2);
axesOpt.btlrMargins(3:4) = axesOpt.btlrMargins(3:4)/axesOpt.figureSize(1);
axesOpt.gapBetweenAxes = axesOpt.gapBetweenAxes./axesOpt.figureSize;

axesHandle = plt.tightSubplot(numOfRows, numOfCols, idx, axesOpt.gapBetweenAxes, axesOpt.btlrMargins(1:2), axesOpt.btlrMargins(3:4));
if axesOpt.reposition; set(gcf, 'position', [screenSize(1:2)+screenSize(3:4)-axesOpt.figureSize-[0 75], axesOpt.figureSize]); end
figureHandle = gcf;
end
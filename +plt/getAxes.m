function [axesHandle, figureHandle] = getAxes(axesOpt)
if ~isfield(axesOpt, 'idx'); axesOpt.idx = 1; end
if ~isfield(axesOpt, 'totalNumOfAxes'); axesOpt.idx = 1; end
if ~isfield(axesOpt, 'figureSize'); axesOpt.figureSize = 400; end
if ~isfield(axesOpt, 'figureHWRatio'); axesOpt.figureHWRatio = 1; end
if ~isfield(axesOpt, 'gapBetweenAxes'); axesOpt.gapBetweenAxes = 25; end
if ~isfield(axesOpt, 'btlrMargins'); axesOpt.btlrMargins = 50*ones(1,4); end
if ~isfield(axesOpt, 'reposition'); axesOpt.reposition = 1; end

screenSize = get(0,'MonitorPositions');
screenSize = screenSize(screenSize(:,1)==min(screenSize(:,1)),:);
screenRatio = round(screenSize(3)/screenSize(4));

if axesOpt.totalNumOfAxes < 4; numOfRows = 1;
else, numOfRows = find(((1:5)*screenRatio.*(1:5))>=axesOpt.totalNumOfAxes,1);
end
numOfCols = ceil(axesOpt.totalNumOfAxes/numOfRows);
axesOpt.figureSize = min([axesOpt.figureSize*numOfCols*axesOpt.figureHWRatio, axesOpt.figureSize*numOfRows], screenSize(3:4));

axesOpt.btlrMargins(1:2) = axesOpt.btlrMargins(1:2)/axesOpt.figureSize(2);
axesOpt.btlrMargins(3:4) = axesOpt.btlrMargins(3:4)/axesOpt.figureSize(1);
axesOpt.gapBetweenAxes = axesOpt.gapBetweenAxes./axesOpt.figureSize;

axesHandle = plt.tightSubplot(numOfRows, numOfCols, axesOpt.idx, axesOpt.gapBetweenAxes, axesOpt.btlrMargins(1:2), axesOpt.btlrMargins(3:4));
if axesOpt.reposition; set(gcf, 'position', [screenSize(1:2)+screenSize(3:4)-axesOpt.figureSize-[0 75], axesOpt.figureSize]); end
figureHandle = gcf;
end
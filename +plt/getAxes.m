function [axesHandle, figureHandle] = getAxes(axIdx, numOfAx, figureSize, figureHWRatio, btlrEdge, axisGap)
if ~exist('axIdx', 'var') || isempty(axIdx); axIdx = 1; end
if ~exist('numOfAx', 'var') || isempty(numOfAx); axIdx = 1; end
if ~exist('figureHWRatio', 'var') || isempty(figureHWRatio); figureHWRatio = 1; end
if ~exist('axisGap', 'var') || isempty(axisGap); axisGap = 25; end
if ~exist('btlrEdge', 'var') || isempty(btlrEdge); btlrEdge = 50*ones(1,4); end
if ~exist('figureSize', 'var') || isempty(figureSize); figureSize = 400; end

screenSize = get(0,'MonitorPositions');
screenSize = screenSize(screenSize(:,1)==min(screenSize(:,1)),:);
screenRatio = round(screenSize(3)/screenSize(4));

if numOfAx < 4; numOfRows = 1;
else, numOfRows = find(((1:5)*screenRatio.*(1:5))>=numOfAx,1);
end
numOfCols = ceil(numOfAx/numOfRows);
figureSize = min([figureSize*numOfCols*figureHWRatio, figureSize*numOfRows], screenSize(3:4));

btlrEdge(1:2) = btlrEdge(1:2)/figureSize(2);
btlrEdge(3:4) = btlrEdge(3:4)/figureSize(1);
axisGap = axisGap./figureSize;

axesHandle = plt.tightSubplot(numOfRows, numOfCols, axIdx, axisGap, btlrEdge(1:2), btlrEdge(3:4));
set(gcf, 'position', [screenSize(1:2)+screenSize(3:4)-figureSize-[0 75], figureSize])
figureHandle = gcf;
end
function viewBoxPlots(obj, plotType, alter)
if ~exist('plotType', 'var'); plotType = 'res'; end
if ~exist('alter', 'var'); alter = 0; figure; else; axesOpt.reposition = 0; end
if isgraphics(alter,'figure'); clf; end

maxGrid = max(cell2mat(cellfun(@(x) [length(x.audValues) length(x.visValues)], obj.blocks, 'uni', 0)), [], 1);
axesOpt.figureHWRatio = maxGrid(2)/(1.3*maxGrid(1));
axesOpt.btlrMargins = [100 80 60 100];
axesOpt.gapBetweenAxes = [100 40];

boxPlot.colorMap = plt.redblue(64);
boxPlot.axisLimits = [0 1];
colorBar.colorLabel = 'Fraction of right turns';
colorBar.colorDirection = 'normal';

if isgraphics(alter,'figure') || ~alter; subjects2Run = 1:length(obj.subjects);
elseif isaxes(alter); disp('NotFunctionalYet');
end

for i  = subjects2Run
    [normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
    boxPlot.subject = obj.subjects{i};
    boxPlot.trialNumber = length(normBlock.responseMade);
    boxPlot.nSessions = obj.blocks{i}.nSessions;
    boxPlot.xyValues = {normBlock.visValues*100; normBlock.audValues};
    boxPlot.xyLabel = {normBlock.audType; 'VisualContrast'};
    switch lower(plotType(1:3))
        case 'res'
            boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade==2, @mean);
            if isempty(obj.figureHandles) || ~any(ismember(obj.figureHandles, gcf)); obj.figureHandles(end+1) = gcf; end
            set(gcf, 'Tag', 'boxRes', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
        case 'gng'
            [~,normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i}, 1, -1);
            normBlock = prc.combineBlocks(normBlock, normBlock.timeOutsBeforeResponse==0);
            set(gcf, 'Tag', 'boxGNG', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
            boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade~=0, @mean);
        case 'las'
            [~,normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i}, 1);
            set(gcf, 'Tag', 'boxLas', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
            boxPlot.plotData = prc.makeGrid(normBlock, normBlock.laserType~=0, @sum);
        case 'num'
            boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade==2, @length);
            set(gcf, 'Tag', 'boxNum', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
            colorBar.colorLabel = 'Relative Num of Trials';
            boxPlot.axisLimits = [0 max(boxPlot.plotData(:))];
        case 'rea'
            boxPlot.plotData = prc.makeGrid(normBlock, round(normBlock.responseTime*1e3), @median, 1);
            boxPlot.axisLimits = [min(boxPlot.plotData(:)) max(boxPlot.plotData(:))];
        case 'tim'
            [normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i}, 0, 5);
            boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade==0, @mean);
            boxPlot.axisLimits = [min(boxPlot.plotData(:)) max(boxPlot.plotData(:))];
    end
    axesOpt.totalNumOfAxes = length(obj.subjects);
    plt.getAxes(axesOpt, i);
    plt.boxPlot(boxPlot);
    colorBar.colorYTick = {'Min'; 'Max'};
end
currentAxisPotision = get(gca, 'position');
figureSize = get(gcf, 'position');

colorBar.handle = colorbar;
set(colorBar.handle, 'Ticks', get(colorBar.handle, 'Limits'), 'TickLabels', colorBar.colorYTick, 'YDir', colorBar.colorDirection);
set(gca, 'position', currentAxisPotision);
colorBar.textHandle = ylabel(colorBar.handle, colorBar.colorLabel);
set(colorBar.textHandle, 'position', [1 mean(get(colorBar.handle, 'Limits')) 0], 'FontSize', 14)
set(colorBar.handle, 'position', [1-75/figureSize(3), 0.2, 30/figureSize(3), 0.6])
end
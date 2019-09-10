function viewBoxPlots(obj, plotType, alter)
if ~exist('plotType', 'var'); plotType = 'res'; end
if ~exist('alter', 'var'); alter = 0; figure; else; axesOpt.reposition = 0; end
if isgraphics(alter,'figure'); clf; end

allBlk = prc.catStructs(obj.blocks);
maxGrid = max(cell2mat(cellfun(@(x) [length(unique(x(:,1))) length(unique(x(:,2)))], allBlk.exp.conditionParametersAV, 'uni', 0)), [], 1);
axesOpt.figureHWRatio = maxGrid(2)/(1.3*maxGrid(1));
axesOpt.btlrMargins = [100 80 60 100];
axesOpt.gapBetweenAxes = [100 40];

boxPlot.colorMap = plt.redblue(64);
boxPlot.axisLimits = [0 1];
colorBar.colorLabel = 'Fraction of right turns';
colorBar.colorDirection = 'normal';

if isgraphics(alter,'figure') || ~alter; subjects2Run = 1:length(obj.blocks);
elseif isaxes(alter); disp('NotFunctionalYet');
end

for i  = subjects2Run
    blk = spatialAnalysis.getBlockType(obj.blocks(i),'norm');
    boxPlot.subject = unique(blk.exp.subject);
    if length(boxPlot.subject) > 1; error('You seem to have a mixed up block file with multiple subjects'); end
    
    boxPlot.trialNumber = blk.tot.trials;
    boxPlot.nExperiments = blk.tot.experiments;
    if boxPlot.nExperiments == 1; boxPlot.extraInf = [blk.exp.expDate{1} ' on ' blk.exp.rigName{1}]; end
    allconditionParametersAV = cell2mat(blk.exp.conditionParametersAV);
    audValues = unique(allconditionParametersAV(:,1));
    visValues = unique(allconditionParametersAV(:,2));
    
    boxPlot.xyValues = {visValues*100; audValues};
    boxPlot.xyLabel = {'AuditoryAzimuth'; 'VisualContrast'};
    switch lower(plotType(1:3))
        case 'res'
            boxPlot.plotData = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @mean);
            if isempty(obj.hand.figure) || ~any(ismember(obj.hand.figure, gcf)); obj.hand.figure(end+1) = gcf; end
            set(gcf, 'Tag', 'boxRes', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
        case 'gng'
            [~,blk] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks(i), 1, -1);
            blk = prc.combineBlocks(blk, blk.timeOutsBeforeResponse==0);
            set(gcf, 'Tag', 'boxGNG', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
            boxPlot.plotData = prc.makeGrid(blk, blk.tri.outcome.responseMade~=0, @mean);
        case 'las'  
            [~,blk] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks(i), 1);
            set(gcf, 'Tag', 'boxLas', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
            boxPlot.plotData = prc.makeGrid(blk, blk.tri.inactivaiton.laserType~=0, @sum);
        case 'num'
            boxPlot.plotData = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @length);
            set(gcf, 'Tag', 'boxNum', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
            colorBar.colorLabel = 'Relative Num of Trials';
            boxPlot.axisLimits = [0 max(boxPlot.plotData(:))];
        case 'rea'
            boxPlot.plotData = prc.makeGrid(blk, round(blk.responseTime*1e3), @median, 1);
            boxPlot.axisLimits = [min(boxPlot.plotData(:)) max(boxPlot.plotData(:))];
        case 'tim'
            [blk] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks(i), 0, 5);
            boxPlot.plotData = prc.makeGrid(blk, blk.tri.outcome.responseMade==0, @mean);
            boxPlot.axisLimits = [min(boxPlot.plotData(:)) max(boxPlot.plotData(:))];
    end
    axesOpt.totalNumOfAxes = length(obj.blocks);
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
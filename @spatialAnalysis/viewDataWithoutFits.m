function viewDataWithoutFits(obj, plotType)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices
%	'rea'--------------------timeToFirstMove vs fration of rightward choices

if ~exist('plotType', 'var'); plotType = 'res'; end
figure;
axesOpt.totalNumOfAxes = length(obj.blks);
axesOpt.btlrMargins = [80 100 80 40];
axesOpt.gapBetweenAxes = [100 60];
axesOpt.numOfRows = ceil(length(obj.blks)/5);
axesOpt.axesSize = [400 450];


allR2 = [];
for i  = 1:length(obj.blks)
    blk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    blk = prc.filtBlock(blk, blk.tri.stim.audAmplitude~=0);
    plotOpt.Marker = '.'; plotOpt.MarkerSize = 20; plotOpt.lineStyle = '-';
    obj.hand.axes = plt.getAxes(axesOpt, i);
    audValues = unique(blk.exp.conditionParametersAV{1}(:,1));
    visValues = unique(blk.exp.conditionParametersAV{1}(:,2));
    
    grds = prc.getGridsFromBlock(blk);
    switch plotType(1:3)
        case 'rea'
            grds.reactionTimeComb(all(isnan(grds.reactionTime),2),:) = [];
            colors =  plt.selectRedBlueColors(grds.audValues(~isinf(grds.audValues(:,1)),1));
            plt.rowsOfGrid(grds.visValues(1,:)*100, grds.reactionTimeComb, colors);
        case 'ken'
            blk = prc.getKennethResopnseFromBlock(blk);
            gridData = prc.makeGrid(blk, round(blk.tri.outcome.timeToKenMove*1e3), @nanmedian, [], 1);
            gridData = nanmean(gridData,3);
            plt.gridSplitByRows(gridData, visValues*100, audValues, plotOpt);
        case 'res'
%             gridData = prc.makeGrid(blk, blk.tri.outcome.responseCalc==2, @mean, 1);
            colors =  plt.selectRedBlueColors(grds.audValues(~isinf(grds.audValues(:,1)),1));
            plt.rowsOfGrid(grds.visValues(1,:)*100, grds.fracRightTurns(1:3,:), colors);
%             plt.gridSplitByRows(gridData, visValues*100, audValues, plotOpt);
            ylim([0 1]);
            xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
            maxContrast = max(abs(blk.tri.stim.visDiff))*100;
            xlim([-maxContrast maxContrast]);
            set(gca, 'xTick', round(((-maxContrast):(maxContrast/4):maxContrast)), 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)));
            yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
        case 'tst'
            blk.tri.outcome.firstMoveDirection = cellfun(@(x) x(1,2), blk.tri.outcome.timeDirAllMoveOnsets);
            gridData = prc.makeGrid(blk, blk.tri.outcome.firstMoveDirection==2, @mean, 1);
            plt.gridSplitByRows(gridData, visValues*100, audValues, plotOpt);
            ylim([0 1]);
            xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
            maxContrast = max(abs(blk.tri.stim.visDiff))*100;
            xlim([-maxContrast maxContrast]);
            set(gca, 'xTick', round(((-maxContrast):(maxContrast/4):maxContrast)), 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)));
            yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
      case 'log'
            gridData = prc.makeGrid(blk, blk.tri.outcome.responseCalc==2, @mean, 1);
            plotOpt.lineStyle = 'none';
            gridData = log((gridData./(1-gridData)));
            plt.gridSplitByRows(gridData, (abs(visValues)*100).^0.7.*sign(visValues), audValues, plotOpt);
            maxContrast = max(abs(blk.tri.stim.visDiff))*100;
            xlim([-maxContrast^0.7 maxContrast^0.7]);
        case 'pro'
            gridData = prc.makeGrid(blk, blk.tri.outcome.responseCalc==2, @mean, 1);
            plotOpt.lineStyle = 'none';
            gridData = norminv(gridData);
            [dataFit] = plt.gridSplitByRows(gridData, visValues*100, audValues, plotOpt);
            maxContrast = max(abs(blk.tri.stim.visDiff))*100;
            xlim([-maxContrast maxContrast]);
            allR2 = [allR2; mean(dataFit.r2)];
    end
    box off;
    title(sprintf('%s: %d Tri', cell2mat(unique(blk.exp.subject)'), blk.tot.trials))
end
figureSize = get(gcf, 'position');
mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
end
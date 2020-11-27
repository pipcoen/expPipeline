function viewConflictChoiceVsRT(obj, plotTag)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices
%	'rea'--------------------timeToFirstMove vs fration of rightward choices

if ~exist('plotTag', 'var'); plotTag = 'kss'; end
xDat = 0.05:0.01:0.5;
[choiceDat] = deal(nan*ones(length(obj.blks), length(xDat)-1));
[histCountsAud, histCountsVis] = deal(nan*ones(length(obj.blks), length(xDat)));
for i  = 1:length(obj.blks)
    blk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    blk = prc.filtBlock(blk, ~isnan(blk.tri.outcome.responseCalc));
    blk = prc.filtBlock(blk, blk.tri.trialType.conflict);
    minTrials = 5;
    
    audSide = (blk.tri.stim.audInitialAzimuth>0)+1;
    responseCalc = blk.tri.outcome.responseCalc;
    time2Move = blk.tri.outcome.reactionTime;
    audChosen = audSide == responseCalc;
    [binCounts, ~, binIdx] = histcounts(time2Move, xDat);
    audChosenForEachRT = arrayfun(@(x) mean(audChosen(binIdx==x)), 1:length(xDat)-1);
    choiceDat(i,binCounts>minTrials) = audChosenForEachRT(binCounts>minTrials);
    histCountsAud(i,:) = ksdensity(time2Move(audChosen), xDat, 'Bandwidth', 0.02);
    histCountsVis(i,:) = ksdensity(time2Move(~audChosen), xDat, 'Bandwidth', 0.02);
end
%%
if strcmpi(plotTag, 'kss')
    cla
    hold on
    plot(xDat, histCountsAud'-histCountsVis', 'color', [0.5 0.5 0.5]);
    plot(xDat, mean(histCountsAud-histCountsVis), 'k', 'linewidth', 3);
    plot(xDat, histCountsAud, 'color',  [0.5, 0, 0.5]);
    plot(xDat, mean(histCountsAud), 'm', 'linewidth', 3);
    plot(xDat, histCountsVis, 'color',  [0 0.5 0.5]);
    plot(xDat, mean(histCountsVis), 'c', 'linewidth', 3);
end
if strcmpi(plotTag, 'ksm')
    cla
    testDat = {histCountsAud; histCountsVis; histCountsAud-histCountsVis};
    meanDataAVC = cell2mat(cellfun(@(x) mean(x), testDat, 'uni', 0));
    stdDataAVC = cell2mat(cellfun(@(x) std(x), testDat, 'uni', 0));
    lowBoundAVC = meanDataAVC-stdDataAVC./sqrt(length(obj.blks));
    upBoundAVC = meanDataAVC+stdDataAVC./sqrt(length(obj.blks));
    
    audPeak = arrayfun(@(x) find(histCountsAud(x,:)==max(histCountsAud(x,:)),1), 1:size(histCountsAud,1));
    visPeak = arrayfun(@(x) find(histCountsVis(x,:)==max(histCountsVis(x,:)),1), 1:size(histCountsAud,1));
    [~, pVal] = ttest(audPeak, visPeak);
    pVal = round(pVal, 2, 'significant');
    title(['P<' num2str(pVal)]);
    
    plotData = cat(3, meanDataAVC, lowBoundAVC, upBoundAVC);
    plotOpt.Marker = 'none';
    colors = [1 0 1; 0 1 1; 0 0 0];
    plt.rowsOfGrid(xDat, plotData(:,:,:), colors, plotOpt);
end
if strcmpi(plotTag, 'fra')
    cla
    xTicks = xDat(2:end)-(diff(xDat(1:2)))/2;
    choiceDatSmth = choiceDat*nan;
    for i = 1:size(choiceDat,1)
        choiceDatSmth(i,:) = smooth(xTicks, choiceDat(i,:), 0.2, 'loess'); 
    end
    micePerPoint = sum(~isnan(choiceDatSmth), 1);
    choicePlot = zeros(1,size(choiceDatSmth,2),3);    
    choicePlot(1,:,2) = nanmean(choiceDatSmth,1);
    choicePlot(1,:,1) = nanmean(choiceDatSmth,1)-(nanstd(choiceDatSmth,1))./sqrt(micePerPoint);
    choicePlot(1,:,3) = nanmean(choiceDatSmth,1)+(nanstd(choiceDatSmth,1))./sqrt(micePerPoint);
    plotOpt.Marker = 'none';
    plt.rowsOfGrid(xTicks, choicePlot, [0 0 0], plotOpt);
end
xlim([0 0.3]);
end
function exampleCellPSTH
%%
s = spatialAnalysis('all', 'm2ephysgood', 1, 0);
kil.loadAtlas;
atlas.tv = tv;
atlas.st = st;
atlas.av = av;
%%
s.blks  = kil.getClusterLoactionsInAllenSpace(s.blks, [], atlas);
s.blks  = prc.filtBlock(s.blks, strcmp(s.blks.clu.parent, 'MOs'));
s.blks  = prc.filtBlock(s.blks,s.blks.tri.trialType.validTrial);
s.blks  = prc.filtBlock(s.blks,s.blks.tri.stim.audAmplitude~=0);
%%
plotCell(s.blks, 'PC045', 166, 95, [0 60]) % Choice cell
%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\4_exampleChoicePSTH', '-pdf', '-painters');
close;

plotCell(s.blks, 'PC045', 166, 97, [0 100]) % AudONOFF cell
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\4_exampleAudOnOffPSTH', '-pdf', '-painters');
close;

plotCell(s.blks, 'PC045', 166, 18, [0 150]) % AudLR cell
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\4_exampleAudLRPSTH', '-pdf', '-painters');
close;

plotCell(s.blks, 'PC045', 42, 6, [0 30]) % VisLR cell
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\4_exampleVisLRPSTH', '-pdf', '-painters');
close;

end

function plotCell(blks, subject, numOfCells, cellNum, yLim)
timeBins = -0.05:0.001:0.2;
42;

axHeight = 250;
axWidth = 250;
nCols = 3;
nRows = 2;
figHeight = nRows*axHeight;
figWidth = nCols*axWidth;

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;

subBlk = prc.filtBlock(blks, strcmp(blks.exp.subject, subject));
cellsPerExp = groupcounts(subBlk.clu.expRef);
subBlk = prc.filtBlock(subBlk, numOfCells==cellsPerExp);
subBlk = prc.filtBlock(subBlk, ~isnan(subBlk.tri.outcome.responseCalc));

for i = 1:3
    if i == 1
        tBlk = prc.filtBlock(subBlk, ~subBlk.tri.trialType.blank & ~subBlk.tri.trialType.auditory);
        eventTimes = tBlk.tri.timeline.visStimPeriodOnOff(:,1);
        rTest = (tBlk.tri.stim.visDiff > 0);
        rControl = (tBlk.tri.outcome.responseCalc == 2)*2;
    elseif i == 2
        tBlk = prc.filtBlock(subBlk, ~subBlk.tri.trialType.blank & ~subBlk.tri.trialType.visual);
        eventTimes = tBlk.tri.timeline.audStimPeriodOnOff(:,1);
        rTest = (tBlk.tri.stim.audDiff > 0);
        rControl = (tBlk.tri.outcome.responseCalc == 2)*2;
    elseif i == 3
        tBlk = prc.filtBlock(subBlk, subBlk.tri.trialType.conflict);
%         tBlk = prc.filtBlock(tBlk, tBlk.tri.stim.visContrast>0.2);
        audStimTimes = tBlk.tri.timeline.audStimPeriodOnOff(:,1);
        visStimTimes = tBlk.tri.timeline.visStimPeriodOnOff(:,1);
        eventTimes = min([audStimTimes visStimTimes],[],2);
        rTest = (tBlk.tri.outcome.responseCalc == 2);
        rControl = (tBlk.tri.stim.audDiff > 0)*2;
    end
    [cluPSTH, timeValues] = kil.getClusterPSTH(tBlk,eventTimes, timeBins);
    filtIdx = num2cell(prc.makeFreqUniform(rControl+rTest,100),1);
    
    rDat = cellfun(@(x) mean(cluPSTH{cellNum+1}(x&rTest,:)), filtIdx, 'uni', 0);
    rDat = smoothdata(mean(cell2mat(rDat')), 'movmean',25);
    
    lDat = cellfun(@(x) mean(cluPSTH{cellNum+1}(x&~rTest,:)), filtIdx, 'uni', 0);
    lDat = smoothdata(mean(cell2mat(lDat')), 'movmean',25);
    
    plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg); cla;
    hold on; box off;
    plot(timeValues, rDat*(1/diff(timeBins(1:2))), '-r', 'linewidth', 2)
    plot(timeValues, lDat*(1/diff(timeBins(1:2))), '-b', 'linewidth', 2)
    ylim(yLim);
    
    plot([0 0], yLim, '--k')
    
    %"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
    %see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
    %     export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_exampleInactModelFits', '-pdf', '-painters');
end
end

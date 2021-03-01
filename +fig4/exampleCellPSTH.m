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
ex
port_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\4_exampleChoicePSTH', '-pdf', '-painters');
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

audStimTimes = subBlk.tri.timeline.audStimPeriodOnOff(:,1);
visStimTimes = subBlk.tri.timeline.visStimPeriodOnOff(:,1);
eventTimes = audStimTimes;%min([audStimTimes visStimTimes],[],2);
for i = 1:3
    if i == 1
        tType = sign(subBlk.tri.stim.visDiff)+2;
        cCol = num2cell('bkr');
    elseif i == 2
        tType = sign(subBlk.tri.stim.audDiff)+2;
        cCol = num2cell('bkr');
    elseif i == 3
        tType = (subBlk.tri.outcome.responseCalc == 2);
        cCol = num2cell('br');
    end
    [cluPSTH, timeValues] = kil.getClusterPSTH(subBlk,eventTimes, timeBins); 
    
    meanDat = arrayfun(@(x) mean(cluPSTH{cellNum+1}(tType==x,:))*(1/diff(timeBins(1:2))), unique(tType), 'uni', 0);
    
    plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg); cla;
    hold on; box off;
    cellfun(@(x,y) plot(timeValues, smoothdata(x,'movmean',25), y),meanDat,cCol');
    ylim(yLim);
    plot([0 0], yLim, '--k')
    
    %"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
    %see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
    %     export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_exampleInactModelFits', '-pdf', '-painters');
end
end

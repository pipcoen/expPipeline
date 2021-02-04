function wheelMovement
%% This function plots the data panels for figure one of the ms
pathInfo.subject = 'PC043';
pathInfo.expDate = '2019-03-22';
s = spatialAnalysis(pathInfo.subject, pathInfo.expDate, 0, 1, 'raw');
pathInfo.expNum = s.blks.exp.expNum{1};
%%
rawBlk = load(prc.pathFinder('backupblock', pathInfo));
rawBlk = rawBlk.block;

% %% Test
% sR = 1000;
% sWin = 51;
% wheelTime = 0:1/sR:rawBlk.inputs.wheelTimes(end);
% wheelDeg = 360*rawBlk.inputs.wheelValues/(4*360);
% wheelDeg = interp1(rawBlk.inputs.wheelTimes, wheelDeg, wheelTime', 'linear', 'extrap')';
% wheelVel = smooth([0 (wheelDeg(2:end)-wheelDeg(1:end-1))*sR],sWin)';
% wheelVel = interp1(wheelTime, wheelVel, wheelTime'-(sWin/2)/sR, 'linear', 'extrap')';
% plot(wheelTime, wheelDeg); hold on;
% plot(wheelTime, wheelDeg+wheelVel)
% 
% blk = s.blks(1); 
% stimTime = blk.tri.timings.stimPeriodStart;
% reactTime = blk.tri.outcome.reactionTime+stimTime;
% gdIdx = strfind((blk.tri.outcome.reactionTime == blk.tri.outcome.timeToFirstMove)', [1 1 1 1 1 1]);
% plotIdx = unique(cell2mat(arrayfun(@(x) x:x+5, gdIdx, 'uni', 0)))';
% plot(reactTime(plotIdx), wheelDeg(round(reactTime(plotIdx)*sR)), '*')
%%

sR = 1000;
sWin = 51;
velThresh = s.blks(1).exp.wheelTicksToDecision{1}*0.2; 

wheelTime = 0:1/sR:rawBlk.inputs.wheelTimes(end);
wheelDeg = 360*rawBlk.inputs.wheelValues/(4*360);
wheelDeg = interp1(rawBlk.inputs.wheelTimes, wheelDeg, wheelTime, 'pchip', 'extrap');

rawVel = [0 (wheelDeg(2:end)-wheelDeg(1:end-1))];
% rawVelShift = interp1(wheelTime, rawVel, wheelTime-(sWin/2)/sR, 'linear', 'extrap');
% dirFilt = ((movmax(rawVelShift>0, sWin)) + (movmax(rawVelShift<0, sWin)))*1;
wheelVelSmth = smooth(rawVel*sR,sWin)';
wheelVel = interp1(wheelTime, wheelVelSmth, wheelTime-floor((sWin/2))/sR, 'linear', 'extrap');
% wheelVelSmth = interp1(wheelTime, wheelVelSmth, wheelTime-(sWin/2)/sR, 'linear', 'extrap');


segTime = 368:1/sR:390;
wheelDegSeg = interp1(wheelTime, wheelDeg, segTime, 'nearest', 'extrap')';
wheelVelSeg = interp1(wheelTime, wheelVel, segTime, 'nearest', 'extrap')';

idx = find(s.blks.tri.timings.stimPeriodStart>segTime(1) & s.blks.tri.timings.stimPeriodStart < segTime(end));


% cla
% reactTimes = s.blks.tri.outcome.reactionTime +s.blks.tri.timings.stimPeriodStart;
% reactTimes = reactTimes(~isnan(reactTimes));
% plot(wheelTime, wheelDeg);
% hold on
% plot(wheelTime, wheelVel+wheelDeg);
% plot(wheelTime, wheelVelSmth+wheelDeg);
% plot(reactTimes, wheelDeg(round(reactTimes*sR)), '*')

reactTimes = s.blks.tri.timings.stimPeriodStart(idx) + s.blks.tri.outcome.reactionTime(idx) - segTime(1);
feedBackTimes = s.blks.tri.timings.stimPeriodStart(idx) + s.blks.tri.outcome.timeToFeedback(idx) - segTime(1);
decThrTimes = s.blks.tri.timings.stimPeriodStart(idx) + s.blks.tri.outcome.timeToResponseThresh(idx) - segTime(1);
stimTimes = s.blks.tri.timings.stimPeriodStart(idx) - segTime(1);

moveDir = s.blks.tri.outcome.responseCalc(idx);
segTime = segTime-segTime(1);
wheelDegSeg = (wheelDegSeg-wheelDegSeg(1))*-1;
wheelVelSeg = wheelVelSeg*-1;
%%
figure;
axHeight = 250;
axWidth = 250;
nCols = 3;
nRows = 2;
figHeight = nRows*axHeight;
figWidth = nCols*axWidth;

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);
%%
axH = plt.tightSubplot(nRows,nCols,1:2,axesGap,botTopMarg,lftRgtMarg); cla
[ax,h1,h2] = plotyy(segTime', wheelDegSeg,segTime', wheelVelSeg);
set(h1, 'Color', 'k')
set(h2, 'Color', 'k')
hold(ax(1), 'on');
hold(ax(2), 'on');
box off;
yL = 150;
velLim = 370;
set(ax(1), 'YLim', [-yL yL], 'XLim', [0 22], 'YTick', [-150 150])
set(ax(2), 'YLim', [-velLim velLim], 'XLim', [0 22], 'YTick', [-velLim velLim])
arrayfun(@(x) patch(ax(1), x+[0 0 0.5 0.5 0], yL*[-1 1 1 -1 -1], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none'), stimTimes);

for i = 1:length(reactTimes)
    moveIdx = (round(reactTimes(i)*sR):round(feedBackTimes(i)*sR))+1;
    if moveDir(i) == 2; lCol = 'r'; else, lCol = 'b'; end
    plot(ax(1), segTime(moveIdx), wheelDegSeg(moveIdx), lCol);
    plot(ax(2), segTime(moveIdx), wheelVelSeg(moveIdx), lCol);
end


plot(ax(2), reactTimes, wheelVelSeg(round(reactTimes*sR)), '.m', 'MarkerSize', 12);
plot(ax(1), decThrTimes, wheelDegSeg(round(decThrTimes*sR)), '.m', 'MarkerSize', 12);
plot(ax(2), xlim, velThresh*[1 1], '--k');
plot(ax(2), xlim, velThresh*[1 1]*-1, '--k');


%%
axH = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); cla; hold on;
blk = s.blks;
blk = prc.filtBlock(blk, ~isnan(blk.tri.outcome.reactionTime));
timeWindow = -0.05:0.001:0.1;

for i=1:blk.tot.trials
    timeRef = timeWindow + blk.tri.outcome.reactionTime(i) +  blk.tri.timings.stimPeriodStart(i);
    wheelPosSeg = interp1(wheelTime, wheelDeg, timeRef, 'nearest', 'extrap')';
    wheelPosPreSeg = interp1(wheelTime, wheelDeg, timeRef-0.01, 'nearest', 'extrap')';
    wheelVelSeg = (wheelPosSeg-wheelPosPreSeg)*100;
%     wheelVelSeg = interp1(wheelTime, wheelVel, timeRef, 'nearest', 'extrap')';
%     rawVelSeg = interp1(wheelTime, rawVel, timeRef, 'nearest', 'extrap')';
%     wheelVelSeg(1:find(timeWindow==0)) = wheelVelSeg(1:find(timeWindow==0)).*(rawVelSeg(1:find(timeWindow==0))~=0);
%     if wheelVelSeg(find(timeWindow==0)+2)==0; keyboard; end
    if blk.tri.outcome.responseCalc(i) == 2; lCol = [1,0,0,0.025]; else, lCol =  [0,0,1,0.025]; end
    plot(timeWindow, wheelVelSeg, 'color', lCol);
end
xlim([timeWindow(1) timeWindow(end)])
ylim([-velLim velLim])
set(gca, 'YTick', [-velLim velLim]);
box off;
plot(xlim, velThresh*[1 1], '--k');
plot(xlim, velThresh*[1 1]*-1, '--k');
%%
axH = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); 
cla; hold on;

opt.Marker = 'none';
load('figSSVMMoveModel_mov', 'svmMod')
cCol = [0 1 1; 1 0 1];

for i = 1:2
    modelPerf = cell2mat(cellfun(@(x) mean(x(:,:,i))', svmMod.modPerf, 'uni', 0))';
    meanData = mean(modelPerf);
    seData = std(modelPerf)./sqrt(size(modelPerf,1));
    plotData = cat(3, meanData, meanData-seData, meanData+seData);
    plt.rowsOfGrid(svmMod.svmTimes, plotData, cCol(i,:), opt);
end
xlim([-0.05 0.1]);
ylim([0.5 1])

axH = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); 
cla; hold on;

load('figSSVMMoveModel_stm', 'svmMod')
for i = 1:2
    modelPerf = cell2mat(cellfun(@(x) mean(x(:,:,i))', svmMod.modPerf, 'uni', 0))';
    meanData = mean(modelPerf);
    seData = std(modelPerf)./sqrt(size(modelPerf,1));;
    plotData = cat(3, meanData, meanData-seData, meanData+seData);
    plt.rowsOfGrid(svmMod.svmTimes, plotData, cCol(i,:), opt);
end
xlim([0.05 0.3]);
ylim([0.5 0.8])
%%
% s = spatialAnalysis('all', 'behaviour', 0, 1, 'raw');
% axH = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); 
% cla; hold on;
% s.viewRightLeftWheelSeparationOverTime;

%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\SupX_wheelMovements', '-pdf', '-painters');
end
function inactivationVsWheel
s = spatialAnalysis('all', 'uniscan', 0, 1, 'raw');
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
timeWindow = -0.01:0.01:0.05;
clear wheelReact wheel;
for i = 1:length(s.blks)
    
    iBlk = s.blks(i);
    mosSites = [0.6 2; 1.8, 2; 0.6, 3; -0.6 2; -1.8, 2; -0.6, 3];
    v1Sites = [1.8 -4; 3,-4; 3,-3; -1.8 -4; -3,-4; -3,-3];
    mosIdx = ismember(iBlk.tri.inactivation.galvoPosition, mosSites, 'rows');
    iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | mosIdx);
    iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime) & iBlk.tri.trialType.validTrial);
    for j = 1:iBlk.tot.trials
        timeRef = timeWindow + iBlk.tri.outcome.reactionTime(j);
        wheelTV = iBlk.tri.raw.wheelTimeValue{j};
        wheelTV(:,2) = 360*wheelTV(:,2)/(4*360)*-1;
        wheelPosSeg = interp1(wheelTV(:,1), wheelTV(:,2), timeRef, 'nearest', 'extrap')';
        wheelReact(j,:) = [0; smooth(diff(wheelPosSeg),1)];
    end
    
    movR = iBlk.tri.outcome.responseCalc == 2;
    lasON = iBlk.tri.inactivation.laserType == 1;
    lasR = iBlk.tri.inactivation.galvoPosition(:,1)>0;
    wheel.movLLasO{i,1} = (wheelReact(~movR & ~lasON,:));
    wheel.movRLasO{i,1} = (wheelReact(movR & ~lasON,:));
    wheel.movLLasL{i,1} = (wheelReact(~movR & lasON & ~lasR,:));
    wheel.movLLasR{i,1} = (wheelReact(~movR & lasON & lasR,:));
    wheel.movRLasL{i,1} = (wheelReact(movR & lasON & ~lasR,:));
    wheel.movRLasR{i,1} = (wheelReact(movR & lasON & lasR,:));
end

figure
hold on;
op2Use = @median;
plot(timeWindow, mean(cell2mat(cellfun(op2Use, wheel.movLLasO, 'uni', 0))), 'r')
plot(timeWindow, mean(cell2mat(cellfun(op2Use, wheel.movRLasO, 'uni', 0))), 'b')
plot(timeWindow, mean(cell2mat(cellfun(op2Use, wheel.movLLasL, 'uni', 0))))
plot(timeWindow,mean(cell2mat(cellfun(op2Use, wheel.movLLasR, 'uni', 0))))
plot(timeWindow, mean(cell2mat(cellfun(op2Use, wheel.movRLasL, 'uni', 0))))
plot(timeWindow, mean(cell2mat(cellfun(op2Use, wheel.movRLasR, 'uni', 0))))
%%


axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('v1');

axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('mos');
%%
% export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\timedInactivations', '-pdf', '-painters');
end
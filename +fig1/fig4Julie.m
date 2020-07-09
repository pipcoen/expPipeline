s = spatialAnalysis('all', 'uniscan', 1, 1);

%%
blks = prc.filtBlock(s.blks, s.blks.tri.outcome.validTrial & s.blks.tri.outcome.responseMade ~=0 & s.blks.tri.stim.audInitialAzimuth==0);
galvoPosition = blks.tri.inactivation.galvoPosition;

filtLMR(:,1) = ismember(galvoPosition, [-3,-2], 'rows');
filtLMR(:,2) = blks.tri.inactivation.laserType == 0;
filtLMR(:,3) = ismember(galvoPosition, [3,-2], 'rows');

colors2Use = 'rkb';
%%
figure;     hold on;
for i = 1:3
    blk = prc.filtBlock(blks, filtLMR(:,i));
    [gridData, gridXY] = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @mean, 1);
    [numData, gridXY] = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @length, 1);
    plotData(i,:) = gridData(2,:);
    trialNum(i,:) = numData(2,:);
    
    
    
    plot(gridXY{2},  plotData(i,:), ['.-' colors2Use(i)]);
end

title(cell2mat(unique(blk.exp.subject)'));
figureSize = get(gcf, 'position');
mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
function scatterPlots(uniBlks)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('uniBlks', 'var'); uniBlks = spatialAnalysis('all', 'uniscan', 0, 1); end
%%
%pre-assign performance and reaction structures with nans
nMice = length(uniBlks.blks);
clear reacN reacL
for i = 1:nMice
    iBlk = prc.filtBlock(uniBlks.blks(i), uniBlks.blks(i).tri.inactivation.galvoPosition(:,2)~=4.5);
    iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
    iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
%     iBlk = prc.filtBlock(iBlk, iBlk.tri.stim.visContrast<0.8);
    
    idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType==1;
    iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
    iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
    iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
    iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
    
    rIdx = (iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0)) & iBlk.tri.inactivation.laserType==0;
    iBlk.tri.stim.audDiff(rIdx) = iBlk.tri.stim.audDiff(rIdx)*-1;
    iBlk.tri.stim.visDiff(rIdx) = iBlk.tri.stim.visDiff(rIdx)*-1;
    iBlk.tri.stim.conditionLabel(rIdx) = -1*iBlk.tri.stim.conditionLabel(rIdx);
          
    tTypeIdx = {iBlk.tri.trialType.auditory;iBlk.tri.trialType.visual; iBlk.tri.trialType.coherent; iBlk.tri.trialType.conflict};
    galvoIdx = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3]};
    
    for j = 1:length(galvoIdx)
        gIdx = ismember(iBlk.tri.inactivation.galvoPosition, galvoIdx{j}, 'rows');
        for k = 1:length(tTypeIdx)
            tBlk = prc.filtBlock(iBlk, (iBlk.tri.inactivation.laserType==0 | gIdx) & tTypeIdx{k});
            normBlk = prc.filtBlock(tBlk, tBlk.tri.inactivation.laserType==0);
            lasBlk = prc.filtBlock(tBlk, tBlk.tri.inactivation.laserType==1);
            lasBlk = prc.filtBlock(lasBlk, lasBlk.tri.stim.visDiff<0 | (lasBlk.tri.stim.visDiff==0 &  lasBlk.tri.stim.audDiff<0));
            
            normGrds = prc.getGridsFromBlock(normBlk);
            lasGrds = prc.getGridsFromBlock(lasBlk);
            
            reacN{i,j,k} = nanmean(normGrds.reactionTimeComb(:))*1000;
            reacL{i,j,k} = nanmean(lasGrds.reactionTimeComb(:))*1000;
        end
    end
end

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

colOrd = 'mykw';
typeOrd = {'Aud'; 'Vis'; 'Coh'; 'Con'};
for i = 1:2
    for j = 1:4
        axH = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
        hold on
        scatter(axH, cell2mat(reacN(:,i,j)), cell2mat(reacL(:,i,j)), 25, colOrd(j), 'filled', 'MarkerEdgeColor', 'k');
        [~, pVal(j,1)] = ttest(cell2mat(reacN(:,i,j)), cell2mat(reacL(:,i,j)));
%         pVal(j,1) = ranksum(cell2mat(reacN(:,i,j)), cell2mat(reacL(:,i,j)));
    end
    xlim([100 420]);
    ylim([100 420]);
    axis square;
    plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
    box off;
    legendCell = arrayfun(@(x,y) [x{1} ': p <' num2str(round(y, 2, 'significant'))], typeOrd,pVal,'uni',0);
    legend(legendCell, 'location', 'none', 'FontSize', 10, 'Position', get(gca, 'Position').*[1, 0.7, 1, 0.3])
    legend('boxoff')
end
%%
% export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_scatterPlots', '-pdf', '-painters');
end
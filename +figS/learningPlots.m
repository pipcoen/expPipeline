function learningPlots
if ~exist('lerBlks', 'var'); lerBlks = spatialAnalysis('all', 'learning', 0, 0, ''); end
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
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 200 figWidth, figHeight]);
%%
axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
numMice = length(lerBlks.blks);
numExp = 35; %hardcoded in prc.keyDates;
expPer = arrayfun(@(x) permute(cell2mat(x.exp.performanceAVM),[3 1 2]), lerBlks.blks, 'uni', 0);
performanceAVM = permute(cell2mat(expPer), [3 2 1]);

stage4Exps = nan*ones(numMice, numExp);
conflictsPresent = arrayfun(@(x) unique(x.tri.expRef(x.tri.trialType.conflict)), lerBlks.blks, 'uni', 0);
for i = 1:numMice; stage4Exps(i, conflictsPresent{i}) = 1; end
expStage = squeeze(sum(~isnan(performanceAVM),1))' + ~isnan(stage4Exps);

stage4Reached = arrayfun(@(x) min(strfind(expStage(x,:)==4, [1 1 1])), 1:numMice)';
stage4Fraction = expStage*0;
for i = 1:numMice; stage4Fraction(i, stage4Reached(i):end) = 1; end

plot(1:length(stage4Fraction), mean(stage4Fraction), 'k', 'linewidth', 2)
xlim([1 30]);
box off;

%%
selectedMouse = 10;
plotOpt.lineStyle = '-';
for i = 1:4
    recType = {'last'; 'last'; 'last'; 'first'};
    expIdx = find(expStage(selectedMouse, :)==i, 2, recType{i});
    disp(expIdx)
    expIdx = ismember(1:length(lerBlks.blks(selectedMouse).exp.subject), expIdx);
    blk = prc.filtBlock(lerBlks.blks(selectedMouse), expIdx, 'exp');
    [~, ~, paramSets] = uniquecell(blk.exp.conditionParametersAV);
    blk = prc.filtBlock(blk, paramSets==mode(paramSets), 'exp');
    blk = prc.filtBlock(blk, ~isnan(blk.tri.outcome.responseCalc));
    grids = prc.getGridsFromBlock(blk);
    plt.tightSubplot(nRows,nCols,i+1,axesGap,botTopMarg,lftRgtMarg); cla;
    plt.rowsOfGrid(grids.visValues(1,:), grids.fracRightTurns, plt.selectRedBlueColors(grids.audValues(:,1)));
    
    xlim([-0.8 0.8]);
    ylim([0 1]);
    xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
    yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
    box off;
    axis square
    ylim([0 1]);
    
    title(['Stage ' num2str(i) ': ' num2str(blk.tot.trials) ' tri, ' num2str(blk.tot.experiments) ' exp'])
end
%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\SupX_mouseLearning', '-pdf', '-painters'); close
end

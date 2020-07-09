function learningPlots

if ~exist('lerBlks', 'var'); lerBlks = spatialAnalysis('all', 'learning', 0, 0, ''); end
figure;
numRow = 1; numCol = 11;
axesHandle = plt.tightSubplot(numRow,numCol,9:11,0.05,[0.15 0.15],[0.05 0.05]);
set(gcf, 'position', get(gcf, 'position').*[2 1 0 0] + [0 0 numCol*100 350]);

numMice = length(lerBlks.blks);
numExp = 35; %hardcoded in prc.keyDates;
expPer = arrayfun(@(x) permute(cell2mat(x.exp.performanceAVM),[3 1 2]), lerBlks.blks, 'uni', 0);
performanceAVM = permute(cell2mat(expPer), [3 2 1]);

stage4Exps = nan*ones(numMice, numExp);
conflictsPresent = arrayfun(@(x) unique(x.tri.expRef(x.tri.trialClass.conflict)), lerBlks.blks, 'uni', 0);
for i = 1:numMice; stage4Exps(i, conflictsPresent{i}) = 1; end
expStage = squeeze(sum(~isnan(performanceAVM),1))' + ~isnan(stage4Exps);

stage4Reached = arrayfun(@(x) min(strfind(expStage(x,:)==4, [1 1 1])), 1:numMice)';
stage4Fraction = expStage*0;
for i = 1:numMice; stage4Fraction(i, stage4Reached(i):end) = 1; end

plot(1:length(stage4Fraction), mean(stage4Fraction), 'k', 'linewidth', 2)
xlim([1 30]);
box off
xlabel('Number of training sessions');
ylabel('Sessions to stage 4');

%%
selectedMouse = 10;
plotOpt.lineStyle = '-';
for i = 1:4
    pltIdx = (i-1)*2+1:(i-1)*2+2;
    recType = {'last'; 'last'; 'last'; 'first'};
    expIdx = find(expStage(selectedMouse, :)==i, 5, recType{i});
    expIdx = ismember(1:length(lerBlks.blks(selectedMouse).exp.subject), expIdx);
    blk = prc.filtBlock(lerBlks.blks(selectedMouse), expIdx, 'exp');
    [~, ~, paramSets] = uniquecell(blk.exp.conditionParametersAV);
    blk = prc.filtBlock(blk, paramSets==mode(paramSets), 'exp');
    plt.tightSubplot(numRow,numCol,pltIdx,0.05,[0.4 0.15],[0.05 0.05]);
    cla;
    grids = prc.getGridsFromBlock(blk);
    plt.fracitonRightChoiceData(grids, plotOpt);
    box off;
    ylim([0 1]);
end
end

% xlabel('Session number');
% ylabel('Performance (%)');
% hold on
% 
% if ~exist('lerBlks', 'var'); lerBlks = spatialAnalysis('all', 'learning', 0, 0, ''); end
% currAx = getPanelAxes(specificPanels, 'b');
% xlabel('Session number');
% ylabel('Performance (%)');
% hold on
% 
% numMice = length(lerBlks.blks);
% numExp = arrayfun(@(x) x.tot.experiments, lerBlks.blks);
% expPer = arrayfun(@(x) permute(cell2mat(x.exp.performanceAVM(1:min(numExp))),[3 1 2]), lerBlks.blks, 'uni', 0);
% performanceAVM = permute(cell2mat(expPer), [3 2 1]);
% 
% performanceAVM > 60;
% 
% mnMul = squeeze(nanmedian(performanceAVM(3,1:15,:),3));
% stMul = squeeze(mad(performanceAVM(3,1:15,:),[],3));
% plt.lineWithPatch(1:15,  [mnMul-stMul; mnMul; mnMul+stMul]);
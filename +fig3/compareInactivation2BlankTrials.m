function scatterPlotsWithLMETest(pType)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
uniBlks = spatialAnalysis('all', 'uniscan', 0, 1);
%%
if ~exist('pType', 'var') || isempty(pType); pType = 'rea'; end

%pre-assign performance and reaction structures with nans
nMice = length(uniBlks.blks);
clear reacN reacL
galvoIdx = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4]};
% tstVal = cell(length(galvoIdx), 2, 2);

for i = 1:nMice
    iBlk = prc.filtBlock(uniBlks.blks(i), uniBlks.blks(i).tri.inactivation.galvoPosition(:,2)~=4.5);
    iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
    
    if strcmpi(pType, 'rea')
        f2Use = 'reactionTime';
        op2use = @median;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
        blankVal(i,1) = nanmedian(iBlk.tri.outcome.reactionTime(iBlk.tri.trialType.blank));
    else
        f2Use = 'responseRecorded';
        op2use = @mean;
        blankVal(i,1) = mean(iBlk.tri.outcome.responseRecorded(iBlk.tri.trialType.blank)==0);
    end
    
    idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType==1;
    iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
    iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
    iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
    iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
    
    for q = 1:2
        tBlk = iBlk;
        if q == 1
            rIdx = (tBlk.tri.stim.audDiff>0 | (tBlk.tri.stim.audDiff==0 & tBlk.tri.stim.visDiff>0)) & tBlk.tri.inactivation.laserType==0;
        else
            rIdx = (tBlk.tri.stim.visDiff>0 | (tBlk.tri.stim.visDiff==0 & tBlk.tri.stim.audDiff>0)) & tBlk.tri.inactivation.laserType==0;
        end
        tBlk.tri.stim.audDiff(rIdx) = tBlk.tri.stim.audDiff(rIdx)*-1;
        tBlk.tri.stim.visDiff(rIdx) = tBlk.tri.stim.visDiff(rIdx)*-1;
        tBlk.tri.stim.conditionLabel(rIdx) = -1*tBlk.tri.stim.conditionLabel(rIdx);
        
        if q == 1
            rIdx = (tBlk.tri.stim.audDiff>0 | (tBlk.tri.stim.audDiff==0 & tBlk.tri.stim.visDiff>0));
        else
            rIdx = (tBlk.tri.stim.visDiff>0 | (tBlk.tri.stim.visDiff==0 & tBlk.tri.stim.audDiff>0));
        end
        tBlk = prc.filtBlock(tBlk, ~rIdx);
        
        tTypeIdx = {tBlk.tri.trialType.coherent; tBlk.tri.trialType.conflict};
        for j = 1:length(galvoIdx)
            gIdx = ismember(tBlk.tri.inactivation.galvoPosition, galvoIdx{j}, 'rows');            
            for k = 1:length(tTypeIdx)
                lasBlk = prc.filtBlock(tBlk, gIdx & tTypeIdx{k});
                lasBlk = prc.filtBlock(lasBlk, lasBlk.tri.inactivation.laserType==1);
                lasBlk.tri.outcome.responseRecorded = lasBlk.tri.outcome.responseRecorded==0;
                
                tstVal(i,j,k,q) = op2use(lasBlk.tri.outcome.(f2Use));
            end
            disp([i q]);
        end
    end
end

tDat = tstVal(:,:,1,2);
[~,~,stats] = anova2([blankVal, tDat]);
multcompare(stats);
figure;
plt.jitter(num2cell(([blankVal, tDat]),1))
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


colOrd = 'kw';
typeOrd = {'Coh'; 'Con'};
siteOrd = {'MOs'; 'V1'; 'A1'};
leftAlign = {'Aud'; 'Vis'};
if strcmpi(pType, 'rea')
    reacN = cellfun(@(x) x*1000, reacN, 'uni', 0);
    reacL = cellfun(@(x) x*1000, reacL, 'uni', 0);
    axLims = {[100 300]; [100 250]; [100 250]};
    sName = 'modifiedReactionTimeScatters';
elseif strcmpi(pType, 'tim')
    axLims = {[0 0.5]; [0 0.5]; [0 0.5]};
    sName = 'modifiedTimeOutScatters';
end
for i = 1:3
    for q = 1:2
        for k = 1:2
            axH = plt.tightSubplot(nRows,nCols,i+(q-1)*3,axesGap,botTopMarg,lftRgtMarg);
            xlim(axLims{i});
            ylim(axLims{i});
            if k == 1; plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2); end
            hold on
            reacRPlt = cell2mat(reacN(:,i,k,q)); reacRPlt(reacRPlt>axLims{i}(2)) = axLims{i}(2);
            reacLPlt = cell2mat(reacL(:,i,k,q)); reacLPlt(reacLPlt>axLims{i}(2)) = axLims{i}(2);
            
            vC = abs(cell2mat(stimC(:,i,k,q)));
            uniV = unique(vC);
            colMod = sqrt(0.8) - sqrt(uniV);
            if k == 1
            arrayfun(@(x,y) scatter(axH, reacRPlt(vC==x), reacLPlt(vC==x), 25, [1,1,1]*y,...
                'filled', 'MarkerEdgeColor', [1,1,1]*y), uniV, colMod); 
            else
                arrayfun(@(x,y) scatter(axH, reacRPlt(vC==x), reacLPlt(vC==x), 25, colOrd(k),...
                    'filled', 'MarkerEdgeColor', [1,1,1]*y), uniV, colMod);
            end
        end
        axis square;
        box off;
        
        pVal = LMEpVal(i,:,q);
        legendCell = arrayfun(@(x,y) [x{1} ': p <' num2str(round(y, 2, 'significant'))], typeOrd,pVal','uni',0);
        legendCell = [legendCell; [siteOrd{i} '-' leftAlign{q} 'Contra']];
        text(axLims{i}(2)*0.6,axLims{i}(2),legendCell, 'fontsize', 10)
    end
end

%%
export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\' sName], '-pdf', '-painters');
end
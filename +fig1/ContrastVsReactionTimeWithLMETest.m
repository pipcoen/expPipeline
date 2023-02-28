function ContrastVsReactionTimeWithLMETest(behBlksOrig)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('behBlksOrig', 'var') || isempty(behBlksOrig); behBlksOrig = spatialAnalysis('all', 'behavior', 0, 1, ''); end
nMice = length(behBlksOrig.blks);
clear reacN LME
stimType = {'Vis'; 'Coh'; 'Con'};
[LME.reacT, LME.stimC, LME.mIdx] = deal(cell(length(stimType),1));

for mouse = 1:nMice
    iBlk = behBlksOrig.blks(mouse);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial & ~iBlk.tri.trialType.blank & ~iBlk.tri.trialType.auditory);
    iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime) & (iBlk.tri.inactivation.laserType==0 | isnan(iBlk.tri.inactivation.laserType)));
    iBlk.tri.datGlobal = iBlk.tri.outcome.reactionTime;

    tBlk = iBlk;
    tTypeIdx = {tBlk.tri.trialType.visual; tBlk.tri.trialType.coherent; tBlk.tri.trialType.conflict};
    for triT = 1:length(tTypeIdx)
        fBlk = tBlk;
        fBlk = prc.filtBlock(fBlk, tTypeIdx{triT});

        blkGrds = prc.getGridsFromBlock(fBlk, 3);
        selGrd = blkGrds.reactionTimeComb;
        globalReacMouseTriT(mouse, triT) = mean(selGrd(~isnan(blkGrds.(f2Use))));

        visVal = blkGrds.visValues(1,:)';
        selGrd = nanmean([selGrd; fliplr(selGrd)],1)';
        tkIdx = ~isnan(selGrd) & visVal > 0;
        LME.reacT{triT} = [LME.reacT{triT}; selGrd(tkIdx)];
        LME.stimC{triT} = [LME.stimC{triT}; visVal(tkIdx)];
        LME.mIdx{triT} = [LME.mIdx{triT}; selGrd(tkIdx)*0+mouse];
        disp(LME.mIdx{triT});
    end
end

%%
LMEtlbs = cellfun(@(w,x,y,z) table(w,x,y, 'VariableNames',{'ReactTime','stimC','Mouse'}),...
    LME.reacT, LME.stimC, LME.mIdx, 'uni', 0);
LMEfits = cellfun(@(x) fitlme(x, 'ReactTime~stimC+(1|Mouse)'), LMEtlbs, 'uni', 0);
pValMouseTriT = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'stimC')), LMEfits);

%%
triType = {'Vis'; 'Coh'; 'Con'};
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

for triT = 1:3
    axH = plt.tightSubplot(nRows,nCols,triT,axesGap,botTopMarg,lftRgtMarg);
    hold on
    for mIdx = 1:nMice
        xDat = LME.stimC{triT}(LME.mIdx{triT}==mIdx);
        yDat = LME.reacT{triT}(LME.mIdx{triT}==mIdx)-globalReacMouseTriT(mIdx);
        plot(xDat, yDat*1000,'color', [0.5 0.5 0.5])
    end

    ylim([-50 40]);
    xlim([0 0.8]);
    box off;

    legendCell = [triType{triT} '  pVal = ' num2str(pValMouseTriT(triT))];
    text(0,max(ylim)*1.1,legendCell, 'fontsize', 10)
end
%     export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\' fName], '-pdf', '-painters');
end

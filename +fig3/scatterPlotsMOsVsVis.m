function scatterPlotsMOsVsVis(opt)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
uniBlks = spatialAnalysis('all', 'uniscan', 0, 1);
%%
if ~exist('opt', 'var'); opt = struct; end
if ~isfield(opt, 'prmType'); opt.prmType = 'rea'; end
if ~isfield(opt, 'pltType'); opt.pltType = 'mean'; end
if ~isfield(opt, 'revNorm'); opt.revNorm = 1; end
if ~isfield(opt, 'siteLoc'); opt.siteLoc = 'contra'; end
if ~isfield(opt, 'offset'); opt.offset = 1; end
%pre-assign performance and reaction structures with nans
nMice = length(uniBlks.blks);
clear reacN reacL
galvoIdx = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4]};
[LME.reacT, LME.stimC, LME.mIdx, LME.lasOn] = deal(cell(length(galvoIdx), 2, 2));

for mouse = 1:nMice
    iBlk = prc.filtBlock(uniBlks.blks(mouse), uniBlks.blks(mouse).tri.inactivation.galvoPosition(:,2)~=4.5);
    iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial & ~iBlk.tri.trialType.blank);
    
    if strcmpi(opt.prmType, 'rea')
        f2Use = 'reactionTimeComb';
        op2use  = @nanmedian;
        iBlk.tri.datGlobal = iBlk.tri.outcome.reactionTime;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime));
    else
        f2Use = 'fracTimeOutComb';
        op2use  = @mean;
        iBlk.tri.datGlobal = iBlk.tri.outcome.responseRecorded==0;
    end
    
    idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType==1;
    iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
    iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
    iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
    iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
    
    for vCon = 1:2
        tBlk = iBlk;
        if opt.revNorm
            if vCon == 1
                rIdx = (tBlk.tri.stim.visDiff>0 | (tBlk.tri.stim.visDiff==0 & tBlk.tri.stim.audDiff>0)) & tBlk.tri.inactivation.laserType==0;
            else
                rIdx = (tBlk.tri.stim.audDiff>0 | (tBlk.tri.stim.audDiff==0 & tBlk.tri.stim.visDiff>0)) & tBlk.tri.inactivation.laserType==0;
            end
            tBlk.tri.stim.audDiff(rIdx) = tBlk.tri.stim.audDiff(rIdx)*-1;
            tBlk.tri.stim.visDiff(rIdx) = tBlk.tri.stim.visDiff(rIdx)*-1;
            tBlk.tri.stim.conditionLabel(rIdx) = -1*tBlk.tri.stim.conditionLabel(rIdx);
        end
        
        if vCon == 1
            rIdx = (tBlk.tri.stim.visDiff>0 | (tBlk.tri.stim.visDiff==0 & tBlk.tri.stim.audDiff>0));
        else
            rIdx = (tBlk.tri.stim.audDiff>0 | (tBlk.tri.stim.audDiff==0 & tBlk.tri.stim.visDiff>0));
        end
        if strcmpi(opt.siteLoc, 'contra')
            tBlk = prc.filtBlock(tBlk, ~rIdx);
        else
            tBlk = prc.filtBlock(tBlk, rIdx);
        end
                tTypeIdx = {tBlk.tri.trialType.coherent; tBlk.tri.trialType.conflict};
        
        for site = 1:length(galvoIdx)
            gIdx = ismember(tBlk.tri.inactivation.galvoPosition,galvoIdx{site}, 'rows');
            for triT = 1:length(tTypeIdx)
                fBlk = tBlk;
                fBlk = prc.filtBlock(fBlk, (fBlk.tri.inactivation.laserType==0 | gIdx) & tTypeIdx{triT});
                normBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==0);
                lasBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==1);
                
                lasGrds = prc.getGridsFromBlock(lasBlk, 3);
                normGrds = prc.getGridsFromBlock(normBlk, 3);
                
%                 globalReacN(mouse,site,triT,vCon) = mean(normGrds.(f2Use)(~isnan(normGrds.(f2Use))&abs(normGrds.visValues)<0.8));
%                 globalReacL(mouse,site,triT,vCon) = mean(lasGrds.(f2Use)(~isnan(lasGrds.(f2Use))&abs(lasGrds.visValues)<0.8));
                
                globalReacN(mouse,site,triT,vCon) = op2use(normBlk.tri.datGlobal);
                globalReacL(mouse,site,triT,vCon) = op2use(lasBlk.tri.datGlobal);
                
                reacN{mouse,site,triT,vCon} = normGrds.(f2Use)(~isnan(normGrds.(f2Use)));
                reacL{mouse,site,triT,vCon} = lasGrds.(f2Use)(~isnan(lasGrds.(f2Use)));
                stimC{mouse,site,triT,vCon} = normGrds.visValues(~isnan(normGrds.(f2Use)));
                
                tkIdx = ~isnan(lasGrds.(f2Use));
                LME.reacT{site,triT,vCon} = [LME.reacT{site,triT,vCon}; lasGrds.(f2Use)(tkIdx)];
                LME.stimC{site,triT,vCon} = [LME.stimC{site,triT,vCon}; lasGrds.visValues(tkIdx)];
                LME.mIdx{site,triT,vCon} = [LME.mIdx{site,triT,vCon}; lasGrds.(f2Use)(tkIdx)*0+mouse];
                LME.lasOn{site,triT,vCon} = [LME.lasOn{site,triT,vCon}; lasGrds.(f2Use)(tkIdx)*0+1];
                
                tkIdx = ~isnan(normGrds.(f2Use));
                LME.reacT{site,triT,vCon} = [LME.reacT{site,triT,vCon}; normGrds.(f2Use)(tkIdx)];
                LME.stimC{site,triT,vCon} = [LME.stimC{site,triT,vCon}; normGrds.visValues(tkIdx)];
                LME.mIdx{site,triT,vCon} = [LME.mIdx{site,triT,vCon}; normGrds.(f2Use)(tkIdx)*0+mouse];
                LME.lasOn{site,triT,vCon} = [LME.lasOn{site,triT,vCon}; normGrds.(f2Use)(tkIdx)*0];
                disp([site triT]);
            end
        end
    end
end

%%
typeOrd = {'Coh'; 'Con'};
if strcmpi(opt.prmType, 'rea')
    if max(vertcat(reacN{:}))<2
        reacN = cellfun(@(x) x*1000, reacN, 'uni', 0);
        reacL = cellfun(@(x) x*1000, reacL, 'uni', 0);
        globalReacN = globalReacN*1000;
        globalReacL = globalReacL*1000;
    end
end

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

axLims = [-50 200];
axH = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
for triT = 1:2
    hold on
    
    DeltaTMOs = globalReacL(:,1,triT,1)- globalReacN(:,1,triT,1);
    DeltaTVis = globalReacL(:,2,triT,1)- globalReacN(:,2,triT,1);
    
    plot([triT-0.25,triT+0.25], [DeltaTMOs DeltaTVis], '-', 'color', [0.5 0.5 0.5]);
    plot([triT-0.25,triT+0.25], [mean(DeltaTMOs), mean(DeltaTVis)], 'k', 'linewidth', 2)
    [~, pVal(triT)] = ttest(DeltaTMOs, DeltaTVis);
end

ylim(axLims)
xlim([0.5 length(typeOrd)+0.5]);
set(gca, 'XTick', 1:length(typeOrd), 'XTickLabel', typeOrd)
box off;
set(gca, 'position', get(gca, 'position').*[1 1 0.6 1])
text(1,max(ylim),{['Coh: ' num2str(pVal(1))]; ['Con: ' num2str(pVal(2))]}, 'fontsize', 10);
%     export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\' fName], '-pdf', '-painters');

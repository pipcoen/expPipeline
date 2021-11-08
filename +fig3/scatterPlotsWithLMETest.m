function scatterPlotsWithLMETest(opt)
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
                
                globalReacN(mouse,site,triT,vCon) = mean(normGrds.(f2Use)(~isnan(normGrds.(f2Use))&abs(normGrds.visValues)<0.8));
                globalReacL(mouse,site,triT,vCon) = mean(lasGrds.(f2Use)(~isnan(lasGrds.(f2Use))&abs(lasGrds.visValues)<0.8));
                
%                 globalReacN(mouse,site,triT,vCon) = op2use(normBlk.tri.datGlobal);
%                 globalReacL(mouse,site,triT,vCon) = op2use(lasBlk.tri.datGlobal);
                
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
LMEtlbs = cellfun(@(w,x,y,z) table(w,z,nominal(abs(x)),y, 'VariableNames',{'ReactTime','LaserOn','vConst','Mouse'}),...
    LME.reacT, LME.stimC, LME.mIdx, LME.lasOn, 'uni', 0);
LMEfits = cellfun(@(x) fitlme(x, 'ReactTime~vConst+LaserOn+(1|Mouse)'), LMEtlbs, 'uni', 0);
LMEpVal = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'LaserOn')), LMEfits);

%%

colOrd = 'kw';
typeOrd = {'Coh'; 'Con'};
siteOrd = {'MOs'; 'V1'; 'A1'};
leftAlign = {'Vis'; 'Aud'};
if strcmpi(opt.prmType, 'rea')
    if max(vertcat(reacN{:}))<2
        reacN = cellfun(@(x) x*1000, reacN, 'uni', 0);
        reacL = cellfun(@(x) x*1000, reacL, 'uni', 0);
        globalReacN = globalReacN*1000;
        globalReacL = globalReacL*1000;
    end
    sName = 'RT_';
elseif strcmpi(opt.prmType, 'tim')
    sName = 'TO_';
end

if strcmpi(opt.pltType, 'mean') || strcmpi(opt.pltType, 'both')
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
    
    if strcmpi(opt.prmType, 'tim')
        axLims = {[0 0.55]; [0 0.55]; [0 0.55]};
        if strcmpi(opt.siteLoc, 'ipsi')
            axLims = {[0 0.55]; [0 0.55]; [0 0.55]};
        end
    else
        axLims = {[-30 175]; [-25 30]; [-25 30]};
    end    
    for site = 1:3
        for vCon = 1:2
            offset = opt.offset * mean(squeeze(globalReacN(:,site,:,vCon)),2);
            for triT = 1:2
                axH = plt.tightSubplot(nRows,nCols,site+(vCon-1)*3,axesGap,botTopMarg,lftRgtMarg);
                hold on
                tDatN = globalReacN(:,site,triT,vCon);
                tDatL = globalReacL(:,site,triT,vCon);
                
                plot([triT-0.25,triT+0.25], [tDatN tDatL]-offset, '-', 'color', [0.5 0.5 0.5]);
                plot([triT-0.25,triT+0.25], [mean(tDatN), mean(tDatL)]-mean(offset), 'k', 'linewidth', 2)
                
            end
            ylim(axLims{site})
            xlim([0.5 length(typeOrd)+0.5]);
            set(gca, 'XTick', 1:length(typeOrd), 'XTickLabel', typeOrd)
            box off;
            set(gca, 'position', get(gca, 'position').*[1 1 0.6 1])
            
            pVal = LMEpVal(site,:,vCon);
            legendCell = arrayfun(@(x,y) [x{1} ': p <' num2str(round(y, 2, 'significant'))], typeOrd,pVal','uni',0);
            legendCell = [legendCell; [siteOrd{site} '-' leftAlign{vCon} 'Contra']];
            text(1,axLims{site}(2),legendCell, 'fontsize', 10)
        end
    end
    fName = [sName 'Mean_nRev' num2str(opt.revNorm) '_' opt.siteLoc];
%     export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\' fName], '-pdf', '-painters');
end

%%

if strcmpi(opt.pltType, 'scat') || strcmpi(opt.pltType, 'both')
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
    
    if strcmpi(opt.prmType, 'tim')
        axLims = {[0 0.5]; [0 0.3]; [0 0.3]};
    else
        axLims = {[100 300]; [100 250]; [100 250]};
    end
    
    for site = 1:3
        for vCon = 1:2
            for triT = 1:2
                
                axH = plt.tightSubplot(nRows,nCols,site+(vCon-1)*3,axesGap,botTopMarg,lftRgtMarg);
                xlim(axLims{site});
                ylim(axLims{site});
                if triT == 1; plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2); end
                hold on
                reacRPlt = cell2mat(reacN(:,site,triT,vCon)); reacRPlt(reacRPlt>axLims{site}(2)) = axLims{site}(2);
                reacLPlt = cell2mat(reacL(:,site,triT,vCon)); reacLPlt(reacLPlt>axLims{site}(2)) = axLims{site}(2);
                
                vC = abs(cell2mat(stimC(:,site,triT,vCon)));
                uniV = unique(vC);
                colMod = sqrt(0.8) - sqrt(uniV);
                if triT == 1
                    arrayfun(@(x,y) scatter(axH, reacRPlt(vC==x), reacLPlt(vC==x), 25, [1,1,1]*y,...
                        'filled', 'MarkerEdgeColor', [1,1,1]*y), uniV, colMod);
                else
                    arrayfun(@(x,y) scatter(axH, reacRPlt(vC==x), reacLPlt(vC==x), 25, colOrd(triT),...
                        'filled', 'MarkerEdgeColor', [1,1,1]*y), uniV, colMod);
                end
            end
            axis square;
            box off;
            
            pVal = LMEpVal(site,:,vCon);
            legendCell = arrayfun(@(x,y) [x{1} ': p <' num2str(round(y, 2, 'significant'))], typeOrd,pVal','uni',0);
            legendCell = [legendCell; [siteOrd{site} '-' leftAlign{vCon} 'Contra']];
            text(axLims{site}(2)*0.6,axLims{site}(2),legendCell, 'fontsize', 10)
        end
    end
    fName = [sName 'Scatter_nRev' num2str(opt.revNorm) '_' opt.siteLoc];
%     export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\' fName], '-pdf', '-painters');
end

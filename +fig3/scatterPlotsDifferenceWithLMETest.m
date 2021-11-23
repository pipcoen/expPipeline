function scatterPlotsDifferenceWithLMETest(opt)
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
[LMESite.reacT, LMESite.stimC, LMESite.mIdx, LMESite.siteIdx] = deal(cell(2, 2));

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
                
                globalReacN(mouse,site,triT,vCon) = mean(normGrds.(f2Use)(~isnan(normGrds.(f2Use))&abs(normGrds.visValues)));
                globalReacL(mouse,site,triT,vCon) = mean(lasGrds.(f2Use)(~isnan(lasGrds.(f2Use))&abs(lasGrds.visValues)));
                
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
                
                tkIdx = ~isnan(lasGrds.(f2Use));
                LMESite.reacT{triT,vCon} = [LMESite.reacT{triT,vCon}; lasGrds.(f2Use)(tkIdx)-normGrds.(f2Use)(tkIdx)];
                LMESite.stimC{triT,vCon} = [LMESite.stimC{triT,vCon}; lasGrds.visValues(tkIdx)];
                LMESite.mIdx{triT,vCon} = [LMESite.mIdx{triT,vCon}; lasGrds.(f2Use)(tkIdx)*0+mouse];
                LMESite.siteIdx{triT,vCon} = [LMESite.siteIdx{triT,vCon}; lasGrds.(f2Use)(tkIdx)*0+site];
                disp([site triT]);
            end
        end
    end
end

%%
LMEtlbs = cellfun(@(w,x,y,z) table(w,z,nominal(abs(x)),y, 'VariableNames',{'ReactTime','LaserOn','vConst','Mouse'}),...
    LME.reacT, LME.stimC, LME.mIdx, LME.lasOn, 'uni', 0);
LMEfits = cellfun(@(x) fitlme(x, 'ReactTime~vConst+LaserOn+(1|Mouse)'), LMEtlbs, 'uni', 0);
LMELaserpVal = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'LaserOn')), LMEfits);

%%
compDo = [1 2; 1,3; 2,3];
for i = 1:3
    tDat = LMESite;
    tDat = structfun(@(x) cellfun(@(y,z) y(ismember(z, compDo(i,:))), x, LMESite.siteIdx, 'uni', 0),tDat, 'uni', 0);
    tDattlbs = cellfun(@(w,x,y,z) table(w,nominal(abs(x)),y,z, 'VariableNames',{'ReactTime','vConst','Mouse','tSite'}),...
        tDat.reacT, tDat.stimC, tDat.mIdx, tDat.siteIdx, 'uni', 0);
    tDatfits = cellfun(@(x) fitlme(x, 'ReactTime~vConst+tSite+(1|Mouse)'), tDattlbs, 'uni', 0);
    LMESitepVal{i,1} = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'tSite')), tDatfits);
end

%%
typeOrd = {'Coh'; 'Con'};
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
        axLims = {[-0.05 0.55]; [-0.05 0.55]; [-0.05 0.55]};
        if strcmpi(opt.siteLoc, 'ipsi')
            axLims = {[-0.05 0.55]; [-0.05 0.55]; [-0.05 0.55]};
        end
    else
        axLims = {[-30 105]; [-30 105]; [-30 105]};
    end    
    for vCon = 1:2
        if vCon == 2
            siteOrd = {'MOs'; 'A1'; 'V1'}; 
            newSiteOrd = [1 3 2];
            newCompOrd = [2 1 3];
        else
            siteOrd = {'MOs'; 'V1'; 'A1'}; 
            newSiteOrd = 1:3;
            newCompOrd = 1:3;
        end
        for triT = 1:2
            axH = plt.tightSubplot(nRows,nCols,triT+(vCon-1)*3,axesGap,botTopMarg,lftRgtMarg);
            hold on
            tDatN = globalReacN(:,:,triT,vCon);
            tDatL = globalReacL(:,:,triT,vCon);
            tDatDiff = tDatL - tDatN;
            tDatDiff(tDatDiff>100) = 100;
            tDatDiff = num2cell(tDatDiff,1);

            nXPnts = length(tDatDiff);
            yDat = cell2mat(arrayfun(@(x) [tDatDiff{x}; mean(tDatDiff{x})], 1:nXPnts, 'uni', 0));
            if vCon == 2 
                yDat(:,1:3) = yDat(:,newSiteOrd);
            end

            xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));
            
            
            set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
            hold on
            for i = 1:nXPnts-1
                cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(1:end-1,i:i+1),2), num2cell(yDat(1:end-1,i:i+1),2));
                cellfun(@(x,y) plot(x,y, 'b','HandleVisibility','off'), num2cell(xDat(end,i:i+1),2), num2cell(yDat(end,i:i+1),2));
            end
            for i = 1:nXPnts
                if LMELaserpVal(newSiteOrd(i),triT,vCon) < 0.05
                    S = num2str(LMELaserpVal(newSiteOrd(i),triT,vCon), '%.1E');
                else
                    S = 'ns';
                end
                    text(xDat(1,i), max(ylim)*1.05, S, ... 
                        'HorizontalAlignment', 'center', 'fontsize', 6);
            end
            
            for i = 1:length(compDo)
                if LMESitepVal{newCompOrd(i)}(triT,vCon) < 0.05
                    S = num2str(LMESitepVal{newCompOrd(i)}(triT,vCon), '%.1E');
                else
                    S = 'ns';
                end
                xPnts = xDat(i, compDo(i,:));
                yLev = max(axLims{i})*(1.05+i*0.05);
                plot(xPnts, yLev*[1 1], 'k')
                text(xPnts(2), yLev, S, ...
                    'HorizontalAlignment', 'left', 'fontsize', 6);
            end
            ylim([axLims{i}(1), yLev]);
            xlim([xDat(1,1)-0.5 xDat(1,end)+0.5]);
            set(gca, 'XTick', xDat(1,:), 'XTickLabel', siteOrd)
            box off;
            
            legendCell = [typeOrd{triT} '-' leftAlign{vCon} 'Contra'];
            text(1,max(ylim)*1.1,legendCell, 'fontsize', 10)
        end
    end
    fName = [sName 'Diff_nRev' num2str(opt.revNorm) '_' opt.siteLoc];
    export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\' fName], '-pdf', '-painters');
end
end

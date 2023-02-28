function scatterPlotsDifferenceBilateralComplete(opt)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
bilBlks = spatialAnalysis('all', 'biscan', 0, 1);
%%
if ~exist('opt', 'var'); opt = struct; end
if ~isfield(opt, 'prmType'); opt.prmType = 'ti2'; end
if ~isfield(opt, 'pltType'); opt.pltType = 2; end
%pre-assign performance and reaction structures with nans
nMice = length(bilBlks.blks);
clear reacN reacL

galvoIdx = {[0.5,2;1.5, 2; 0.5, 3];[1.5 -4; 2.5,-4; 2.5,-3];[1.5,-2; 2.5,-2; 3.5,-2]};
[LME.reacT, LME.stimC, LME.mIdx, LME.lasOn] = deal(cell(length(galvoIdx), 4));
[LMESite.reacT, LMESite.stimC, LMESite.mIdx, LMESite.siteIdx] = deal(cell(4, 1));

for mouse = 1:nMice
    iBlk = prc.filtBlock(bilBlks.blks(mouse), bilBlks.blks(mouse).tri.inactivation.galvoPosition(:,2)~=4.5);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial & ~iBlk.tri.trialType.blank);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.stim.visContrast < 0.07);
    iBlk.tri.stim.audDiff(isinf(iBlk.tri.stim.audDiff)) = 0;

    rtLimit = 1.5;
    iBlk.tri.outcome.responseCalc(iBlk.tri.outcome.reactionTime>rtLimit) = nan;
    iBlk.tri.outcome.responseRecorded(iBlk.tri.outcome.reactionTime>rtLimit) = 0;
    iBlk.tri.outcome.reactionTime(iBlk.tri.outcome.reactionTime>rtLimit) = nan;

    rIdx = iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0);
    iBlk.tri.outcome.responseCalc(rIdx) = (iBlk.tri.outcome.responseCalc(rIdx)*-1+3).*(iBlk.tri.outcome.responseCalc(rIdx)>0);
    iBlk.tri.stim.audInitialAzimuth(rIdx) = iBlk.tri.stim.audInitialAzimuth(rIdx)*-1;
    iBlk.tri.stim.visInitialAzimuth(rIdx) = iBlk.tri.stim.visInitialAzimuth(rIdx)*-1;
    iBlk.tri.stim.visInitialAzimuth(isinf(iBlk.tri.stim.visInitialAzimuth)) = inf;
    iBlk.tri.stim.conditionLabel(rIdx) = iBlk.tri.stim.conditionLabel(rIdx)*-1;

    
    if strcmpi(opt.prmType, 'res')
        f2Use = 'fracRightTurnsComb';
        iBlk.tri.datGlobal = iBlk.tri.outcome.responseCalc;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
    elseif strcmpi(opt.prmType, 'rea')
        f2Use = 'reactionTimeComb';
        iBlk.tri.datGlobal = iBlk.tri.outcome.reactionTime;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime));
    elseif strcmpi(opt.prmType, 'tim')
        f2Use = 'fracTimeOutComb';
        iBlk.tri.datGlobal = iBlk.tri.outcome.responseRecorded==0;
    elseif strcmpi(opt.prmType, 'ti2')
        f2Use = 'fracLongResponses';
        iBlk.tri.datGlobal = iBlk.tri.outcome.reactionTime;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime));
    end

    tTypeIdx = {iBlk.tri.trialType.visual; iBlk.tri.trialType.auditory; iBlk.tri.trialType.coherent; iBlk.tri.trialType.conflict};
    for site = 1:length(galvoIdx)
        gPosAbs = [abs(iBlk.tri.inactivation.galvoPosition(:,1)) iBlk.tri.inactivation.galvoPosition(:,2)];
        gIdx = ismember(gPosAbs,galvoIdx{site}, 'rows');
        for triT = 1:length(tTypeIdx)
            fBlk = iBlk;
            fBlk = prc.filtBlock(fBlk, (fBlk.tri.inactivation.laserType==0 | gIdx) & tTypeIdx{triT});

            normBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==0);
            lasBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==2);

            lasGrds = prc.getGridsFromBlock(lasBlk, 3);
            normGrds = prc.getGridsFromBlock(normBlk, 3);

            lasGrds.visValues(lasGrds.visValues==-0.06) = -0.05;
            lasGrds.visValues(lasGrds.visValues==0.06) = 0.05;
            normGrds.visValues(normGrds.visValues==-0.06) = -0.05;
            normGrds.visValues(normGrds.visValues==0.06) = 0.05;

            globalReacN(mouse,site,triT) = mean(normGrds.(f2Use)(~isnan(normGrds.(f2Use))));
            globalReacL(mouse,site,triT) = mean(lasGrds.(f2Use)(~isnan(lasGrds.(f2Use))));

            reacN{mouse,site,triT} = normGrds.(f2Use)(~isnan(normGrds.(f2Use)));
            reacL{mouse,site,triT} = lasGrds.(f2Use)(~isnan(lasGrds.(f2Use)));
            stimC{mouse,site,triT} = normGrds.visValues(~isnan(normGrds.(f2Use)));

            tkIdx = ~isnan(lasGrds.(f2Use));
            lasGrds.visValues(lasGrds.visValues==0) = [-10; 0; 10];
            LME.reacT{site,triT} = [LME.reacT{site,triT}; lasGrds.(f2Use)(tkIdx)];
            LME.stimC{site,triT} = [LME.stimC{site,triT}; lasGrds.visValues(tkIdx)];
            LME.mIdx{site,triT} = [LME.mIdx{site,triT}; lasGrds.(f2Use)(tkIdx)*0+mouse];
            LME.lasOn{site,triT} = [LME.lasOn{site,triT}; lasGrds.(f2Use)(tkIdx)*0+1];

            tkIdx = ~isnan(normGrds.(f2Use));
            normGrds.visValues(normGrds.visValues==0) = [-10; 0; 10];
            LME.reacT{site,triT} = [LME.reacT{site,triT}; normGrds.(f2Use)(tkIdx)];
            LME.stimC{site,triT} = [LME.stimC{site,triT}; normGrds.visValues(tkIdx)];
            LME.mIdx{site,triT} = [LME.mIdx{site,triT}; normGrds.(f2Use)(tkIdx)*0+mouse];
            LME.lasOn{site,triT} = [LME.lasOn{site,triT}; normGrds.(f2Use)(tkIdx)*0];
            disp([site triT]);

            tkIdx = ~isnan(lasGrds.(f2Use));
            LMESite.reacT{triT} = [LMESite.reacT{triT}; lasGrds.(f2Use)(tkIdx)-normGrds.(f2Use)(tkIdx)];
            LMESite.stimC{triT} = [LMESite.stimC{triT}; lasGrds.visValues(tkIdx)];
            LMESite.mIdx{triT} = [LMESite.mIdx{triT}; lasGrds.(f2Use)(tkIdx)*0+mouse];
            LMESite.siteIdx{triT} = [LMESite.siteIdx{triT}; lasGrds.(f2Use)(tkIdx)*0+site];
            disp([site triT]);
        end
    end
end

%%
LMEtlbs = cellfun(@(w,x,y,z) table(w,z,nominal(x),y, 'VariableNames',{'ReactTime','LaserOn','stimC','Mouse'}),...
    LME.reacT, LME.stimC, LME.mIdx, LME.lasOn, 'uni', 0);
LMEfits = cellfun(@(x) fitlme(x, 'ReactTime~LaserOn+stimC+(1|Mouse)'), LMEtlbs, 'uni', 0);
pValLas_Site_triT = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'LaserOn')), LMEfits);


%%
compDo = [1 2; 1,3; 2,3];
for i = 1:3
    tDat = LMESite;
    tDat = structfun(@(x) cellfun(@(y,z) y(ismember(z, compDo(i,:))), x, LMESite.siteIdx, 'uni', 0),tDat, 'uni', 0);
    tDattlbs = cellfun(@(w,x,y,z) table(w,nominal(abs(x)),y,z, 'VariableNames',{'ReactTime','stimC','Mouse','tSite'}),...
        tDat.reacT, tDat.stimC, tDat.mIdx, tDat.siteIdx, 'uni', 0);
    tDatfits = cellfun(@(x) fitlme(x, 'ReactTime~stimC+tSite+(1|Mouse)'), tDattlbs, 'uni', 0);
    pValSite_triT{i,1} = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'tSite')), tDatfits);
end

%%
typeOrd = {'Vis'; 'Aud'; 'Coh'; 'Con'};
leftAlign = {'Contra'; 'Ipsi'};
siteOrd = {'MOs'; 'Vis'; 'PPC'};
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
elseif strcmpi(opt.prmType, 'res')
    sName = 'FR_';
elseif strcmpi(opt.prmType, 'ti2')
    sName = 'LR_';
end
axLims = [-35 155];
if strcmpi(opt.prmType, 'tim'); axLims = [-0.08 0.25]; end
if strcmpi(opt.prmType, 'ti2'); axLims = [-0.05 0.25]; end
if strcmpi(opt.prmType, 'res'); axLims = [-0.15 0.40]; end

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

if opt.pltType == 1
    for targSite = 1:3
        axH = plt.tightSubplot(nRows,nCols,targSite,axesGap,botTopMarg,lftRgtMarg);
        hold on
        tDatN = globalReacN(:,targSite,:);
        tDatL = globalReacL(:,targSite,:);
        tDatDiff = tDatL - tDatN;
        tDatDiff(tDatDiff>100) = 100;
        tDatDiff = num2cell(tDatDiff,1);

        nXPnts = length(tDatDiff);
        yDat = cell2mat(arrayfun(@(x) [tDatDiff{x}; mean(tDatDiff{x})], 1:nXPnts, 'uni', 0));
        xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));


        set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
        hold on
        for i = 1:nXPnts-1
            cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(1:end-1,i:i+1),2), num2cell(yDat(1:end-1,i:i+1),2));
            cellfun(@(x,y) plot(x,y, 'b','HandleVisibility','off'), num2cell(xDat(end,i:i+1),2), num2cell(yDat(end,i:i+1),2));
        end
        for i = 1:nXPnts
            if pValLas_Site_triT(targSite,i) < 0.05
                S = num2str(pValLas_Site_triT(targSite,i), '%.1E');
            else
                S = 'ns';
            end
            text(xDat(1,i), max(ylim)*1.05, S, ...
                'HorizontalAlignment', 'center', 'fontsize', 8);
        end

        ylim(axLims);
        xlim([xDat(1,1)-0.5 xDat(1,end)+0.5]);
        set(gca, 'XTick', xDat(1,:), 'XTickLabel', typeOrd)
        box off;

        legendCell = [siteOrd{targSite} 'Trials-leftStim'];
        text(1,max(ylim)*1.1,legendCell, 'fontsize', 10)
    end
    fName = [sName 'Alt_InactiveDelta_Bilateral'];
%     export_fig(['D:\OneDrive - University College London\Papers\Coen_2021\NeuronRevision\NewFigParts' fName], '-pdf', '-painters');
end


if opt.pltType == 2
    for triType = 1:length(typeOrd)
        axH = plt.tightSubplot(nRows,nCols,triType,axesGap,botTopMarg,lftRgtMarg);
        hold on
        tDatN = globalReacN(:,:,triType);
        tDatL = globalReacL(:,:,triType);
        tDatDiff = tDatL - tDatN;
        tDatDiff(tDatDiff>150) = 150;
        tDatDiff = num2cell(tDatDiff,1);

        nXPnts = length(tDatDiff);
        yDat = cell2mat(arrayfun(@(x) [tDatDiff{x}; mean(tDatDiff{x})], 1:nXPnts, 'uni', 0));
        xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));


        set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
        hold on
        for i = 1:nXPnts-1
            cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(1:end-1,i:i+1),2), num2cell(yDat(1:end-1,i:i+1),2));
            cellfun(@(x,y) plot(x,y, 'b','HandleVisibility','off'), num2cell(xDat(end,i:i+1),2), num2cell(yDat(end,i:i+1),2));
        end
        for i = 1:nXPnts
            if pValLas_Site_triT(i,triType) < 0.05
                S = num2str(pValLas_Site_triT(i,triType), '%.1E');
            else
                S = 'ns';
            end
            text(xDat(1,i), max(ylim)*1.05, S, ...
                'HorizontalAlignment', 'center', 'fontsize', 8);
        end

        for i = 1:length(compDo)
            if pValSite_triT{i}(triType) < 0.05
                S = num2str(pValSite_triT{i}(triType), '%.1E');
            else
                S = 'ns';
            end
            xPnts = xDat(i, compDo(i,:));
            yLev = max(axLims)*(1.05+i*0.05);
            plot(xPnts, yLev*[1 1], 'k')
            text(xPnts(2), yLev, S, ...
                'HorizontalAlignment', 'left', 'fontsize', 8);
        end
        ylim([axLims(1) yLev]);
        xlim([xDat(1,1)-0.5 xDat(1,end)+0.5]);
        set(gca, 'XTick', xDat(1,:), 'XTickLabel', siteOrd)
        box off;

        legendCell = [typeOrd{triType} ' Trials-leftStim'];
        text(1,max(ylim)*1.1,legendCell, 'fontsize', 10)
    end
    fName = [sName 'Alt_InactiveDelta_Bilateral'];
%     export_fig(['D:\OneDrive - University College London\Papers\Coen_2021\NeuronRevision\NewFigParts\' fName], '-pdf', '-painters');
end

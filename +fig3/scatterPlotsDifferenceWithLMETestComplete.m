function scatterPlotsDifferenceWithLMETestComplete(opt)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
uniBlks = spatialAnalysis('all', 'uniscan', 0, 1);
%%
if ~exist('opt', 'var'); opt = struct; end
if ~isfield(opt, 'prmType'); opt.prmType = 'ti2'; end
if ~isfield(opt, 'pltType'); opt.pltType = 2; end
%pre-assign performance and reaction structures with nans
nMice = length(uniBlks.blks);
clear reacN reacL
galvoIdx = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4]};
[LME.reacT, LME.stimC, LME.mIdx, LME.lasOn] = deal(cell(length(galvoIdx), 4, 2));
[LMESite.reacT, LMESite.stimC, LMESite.mIdx, LMESite.siteIdx] = deal(cell(4, 2));

for mouse = 1:nMice
    iBlk = prc.filtBlock(uniBlks.blks(mouse), uniBlks.blks(mouse).tri.inactivation.galvoPosition(:,2)~=4.5);
    iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial & ~iBlk.tri.trialType.blank);
    
    if strcmpi(opt.prmType, 'rea')
        f2Use = 'reactionTimeComb';
        iBlk.tri.datGlobal = iBlk.tri.outcome.reactionTime;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime));
    elseif strcmpi(opt.prmType, 'tim')
        f2Use = 'fracTimeOutComb';
        iBlk.tri.datGlobal = iBlk.tri.outcome.responseRecorded==0;
    elseif strcmpi(opt.prmType, 'ti2')
        f2Use = 'fracLongResponses';
%         iBlk.tri.outcome.reactionTime(iBlk.tri.outcome.responseRecorded == 0) = 1.5;
        iBlk.tri.datGlobal = iBlk.tri.outcome.reactionTime;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime));
    end
    
    for Cn_Ip = 1:2
        tBlk = iBlk;
        tTypeIdx = {tBlk.tri.trialType.visual; tBlk.tri.trialType.auditory; tBlk.tri.trialType.coherent; tBlk.tri.trialType.conflict};
        for site = 1:length(galvoIdx)
            gPosAbs = [abs(tBlk.tri.inactivation.galvoPosition(:,1)) tBlk.tri.inactivation.galvoPosition(:,2)];
            gIdx = ismember(gPosAbs,galvoIdx{site}, 'rows');
            for triT = 1:length(tTypeIdx)
                fBlk = tBlk;
                fBlk = prc.filtBlock(fBlk, (fBlk.tri.inactivation.laserType==0 | gIdx) & tTypeIdx{triT});
                
                normBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==0);
                lasBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==1);
                
                gPos = lasBlk.tri.inactivation.galvoPosition(:,1);
                aDiff = lasBlk.tri.stim.audDiff;
                vDiff = lasBlk.tri.stim.visDiff;
                if Cn_Ip == 1
                    cIdx = sign(vDiff).*sign(gPos) < 0 | (vDiff==0 & (sign(aDiff).*sign(gPos)<0));
                elseif Cn_Ip == 2
                    cIdx = sign(vDiff).*sign(gPos) > 0 | (vDiff==0 & (sign(aDiff).*sign(gPos)>0));
                end
                lasBlk = prc.filtBlock(lasBlk, cIdx);
                
                lasGrds = prc.getGridsFromBlock(lasBlk, 3);
                normGrds = prc.getGridsFromBlock(normBlk, 3);
                
                globalReacN(mouse,site,triT,Cn_Ip) = mean(normGrds.(f2Use)(~isnan(normGrds.(f2Use))));
                globalReacL(mouse,site,triT,Cn_Ip) = mean(lasGrds.(f2Use)(~isnan(lasGrds.(f2Use))));
                
                reacN{mouse,site,triT,Cn_Ip} = normGrds.(f2Use)(~isnan(normGrds.(f2Use)));
                reacL{mouse,site,triT,Cn_Ip} = lasGrds.(f2Use)(~isnan(lasGrds.(f2Use)));
                stimC{mouse,site,triT,Cn_Ip} = normGrds.visValues(~isnan(normGrds.(f2Use)));
                
                tkIdx = ~isnan(lasGrds.(f2Use));
                lasGrds.visValues(lasGrds.visValues==0) = [-10; 0; 10];
                LME.reacT{site,triT,Cn_Ip} = [LME.reacT{site,triT,Cn_Ip}; lasGrds.(f2Use)(tkIdx)];
                LME.stimC{site,triT,Cn_Ip} = [LME.stimC{site,triT,Cn_Ip}; lasGrds.visValues(tkIdx)];
                LME.mIdx{site,triT,Cn_Ip} = [LME.mIdx{site,triT,Cn_Ip}; lasGrds.(f2Use)(tkIdx)*0+mouse];
                LME.lasOn{site,triT,Cn_Ip} = [LME.lasOn{site,triT,Cn_Ip}; lasGrds.(f2Use)(tkIdx)*0+1];
                
                tkIdx = ~isnan(normGrds.(f2Use));
                normGrds.visValues(normGrds.visValues==0) = [-10; 0; 10];
                LME.reacT{site,triT,Cn_Ip} = [LME.reacT{site,triT,Cn_Ip}; normGrds.(f2Use)(tkIdx)];
                LME.stimC{site,triT,Cn_Ip} = [LME.stimC{site,triT,Cn_Ip}; normGrds.visValues(tkIdx)];
                LME.mIdx{site,triT,Cn_Ip} = [LME.mIdx{site,triT,Cn_Ip}; normGrds.(f2Use)(tkIdx)*0+mouse];
                LME.lasOn{site,triT,Cn_Ip} = [LME.lasOn{site,triT,Cn_Ip}; normGrds.(f2Use)(tkIdx)*0];
                disp([site triT]);
                
                tkIdx = ~isnan(lasGrds.(f2Use));
                LMESite.reacT{triT,Cn_Ip} = [LMESite.reacT{triT,Cn_Ip}; lasGrds.(f2Use)(tkIdx)-normGrds.(f2Use)(tkIdx)];
                LMESite.stimC{triT,Cn_Ip} = [LMESite.stimC{triT,Cn_Ip}; lasGrds.visValues(tkIdx)];
                LMESite.mIdx{triT,Cn_Ip} = [LMESite.mIdx{triT,Cn_Ip}; lasGrds.(f2Use)(tkIdx)*0+mouse];
                LMESite.siteIdx{triT,Cn_Ip} = [LMESite.siteIdx{triT,Cn_Ip}; lasGrds.(f2Use)(tkIdx)*0+site];
                disp([site triT]);
            end
        end
    end
end

%%
LMEtlbs = cellfun(@(w,x,y,z) table(w,z,nominal(x),y, 'VariableNames',{'ReactTime','LaserOn','stimC','Mouse'}),...
    LME.reacT, LME.stimC, LME.mIdx, LME.lasOn, 'uni', 0);
LMEfits = cellfun(@(x) fitlme(x, 'ReactTime~stimC+LaserOn+(1|Mouse)'), LMEtlbs, 'uni', 0);
pValLas_Site_triT_Cn_Ip = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'LaserOn')), LMEfits);

%%
% tDat = struct;
% tDat.reacT = [];
% tDat.Cn_Ip = [];
% tDat.mIdx = [];
% tDat.stimC = [];
% for i = 1:3
%     for j = 1:2
%     sIdx = (LMESite.siteIdx{i,j}==1);
%     tDat.reacT = [tDat.reacT; LMESite.reacT{i,j}(sIdx)];
%     tDat.mIdx = [tDat.mIdx; LMESite.mIdx{i,j}(sIdx)];
%     tDat.stimC = [tDat.stimC; LMESite.stimC{i,j}(sIdx)];
%     tDat.Cn_Ip = [tDat.Cn_Ip; LMESite.stimC{i,j}(sIdx)*0+j];
%     end
% end
% LMEtlbs = table(tDat.reacT,tDat.Cn_Ip,nominal(tDat.stimC),tDat.mIdx, 'VariableNames',{'ReactTime','Cn_Ip','stimC','Mouse'});
% LMEfits = fitlme(LMEtlbs, 'ReactTime~stimC+Cn_Ip+(1|Mouse)');
% pValLas_Site_triT_Cn_Ip = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'LaserOn')), LMEfits);

%%
compDo = [1 2; 1,3; 2,3];
for i = 1:3
    tDat = LMESite;
    tDat = structfun(@(x) cellfun(@(y,z) y(ismember(z, compDo(i,:))), x, LMESite.siteIdx, 'uni', 0),tDat, 'uni', 0);
    tDattlbs = cellfun(@(w,x,y,z) table(w,nominal(abs(x)),y,z, 'VariableNames',{'ReactTime','stimC','Mouse','tSite'}),...
        tDat.reacT, tDat.stimC, tDat.mIdx, tDat.siteIdx, 'uni', 0);
    tDatfits = cellfun(@(x) fitlme(x, 'ReactTime~stimC+tSite+(1|Mouse)'), tDattlbs, 'uni', 0);
    pValSite_triT_Cn_Ip{i,1} = cellfun(@(x) x.Coefficients.pValue(contains(x.Coefficients.Name, 'tSite')), tDatfits);
end

%%
typeOrd = {'Vis'; 'Aud'; 'Coh'; 'Con'};
leftAlign = {'Contra'; 'Ipsi'};
siteOrd = {'MOs'; 'Vis'; 'Aud'};
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
axLims = [-35 105];
if strcmpi(opt.prmType, 'tim'); axLims = [-0.08 0.45]; end
if strcmpi(opt.prmType, 'ti2'); axLims = [-0.2 0.45]; end

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
    for Cn_Ip = 1:2
        for targSite = 1:3
            axH = plt.tightSubplot(nRows,nCols,targSite+(Cn_Ip-1)*3,axesGap,botTopMarg,lftRgtMarg);
            hold on
            tDatN = globalReacN(:,targSite,:,Cn_Ip);
            tDatL = globalReacL(:,targSite,:,Cn_Ip);
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
                if pValLas_Site_triT_Cn_Ip(targSite,i,Cn_Ip) < 0.05
                    S = num2str(pValLas_Site_triT_Cn_Ip(targSite,i,Cn_Ip), '%.1E');
                else
                    S = 'ns';
                end
                text(xDat(1,i), max(ylim)*1.05, S, ...
                    'HorizontalAlignment', 'center', 'fontsize', 6);
            end
      
            ylim(axLims);
            xlim([xDat(1,1)-0.5 xDat(1,end)+0.5]);
            set(gca, 'XTick', xDat(1,:), 'XTickLabel', typeOrd)
            box off;
            
            legendCell = [siteOrd{targSite} ' Inativation-' leftAlign{Cn_Ip}];
            text(1,max(ylim)*1.1,legendCell, 'fontsize', 10)
        end
    end
%     fName = [sName 'InactiveDelta'];
%     export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\' fName], '-pdf', '-painters');
end
%%
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);
if opt.pltType == 2
    for Cn_Ip = 1:2
        figure;
        for triType = 1:length(typeOrd)
            axH = plt.tightSubplot(nRows,nCols,triType,axesGap,botTopMarg,lftRgtMarg);
            hold on
            tDatN = globalReacN(:,:,triType,Cn_Ip);
            tDatL = globalReacL(:,:,triType,Cn_Ip);
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
                if pValLas_Site_triT_Cn_Ip(i,triType,Cn_Ip) < 0.05
                    S = num2str(pValLas_Site_triT_Cn_Ip(i,triType,Cn_Ip), '%.1E');
                else
                    S = 'ns';
                end
                text(xDat(1,i), max(ylim)*1.05, S, ...
                    'HorizontalAlignment', 'center', 'fontsize', 6);
            end
            
            for i = 1:length(compDo)
                if pValSite_triT_Cn_Ip{i}(triType,Cn_Ip) < 0.05
                    S = num2str(pValSite_triT_Cn_Ip{i}(triType,Cn_Ip), '%.1E');
                else
                    S = 'ns';
                end
                xPnts = xDat(i, compDo(i,:));
                yLev = max(axLims)*(1.05+i*0.05);
                plot(xPnts, yLev*[1 1], 'k')
                text(xPnts(2), yLev, S, ...
                    'HorizontalAlignment', 'left', 'fontsize', 6);
            end
            ylim([axLims(1) yLev]);
            xlim([xDat(1,1)-0.5 xDat(1,end)+0.5]);
            set(gca, 'XTick', xDat(1,:), 'XTickLabel', siteOrd)
            box off;
            
            legendCell = [typeOrd{triType} ' Trials-' leftAlign{Cn_Ip}];
            text(1,max(ylim)*1.1,legendCell, 'fontsize', 10)
        end
%         fName = [sName 'Alt_InactiveDelta_' leftAlign{Cn_Ip}];
        export_fig(['C:\Users\Pip\OneDrive - University College London\Papers\Coen_2021\NeuronRevision\Round2\NewFigures\Rev_SlowResponses_' leftAlign{Cn_Ip}], '-pdf', '-painters');
    end
end

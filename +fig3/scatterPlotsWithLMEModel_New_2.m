function scatterPlotsWithLMEModel_New_2(uniBlks)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('uniBlks', 'var'); uniBlks = spatialAnalysis('all', 'uniscan', 0, 1); end

%pre-assign performance and reaction structures with nans

clear tblData;
nMice = length(uniBlks.blks);
for i = 1:nMice
    iBlk = prc.filtBlock(uniBlks.blks(i), uniBlks.blks(i).tri.inactivation.galvoPosition(:,2)~=4.5);
    iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
    iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
    
    idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType==1;
    iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
    iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
    iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
    iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
    
    rIdx = (iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0)) & iBlk.tri.inactivation.laserType==0;
    iBlk.tri.stim.audDiff(rIdx) = iBlk.tri.stim.audDiff(rIdx)*-1;
    iBlk.tri.stim.visDiff(rIdx) = iBlk.tri.stim.visDiff(rIdx)*-1;
    iBlk.tri.stim.conditionLabel(rIdx) = -1*iBlk.tri.stim.conditionLabel(rIdx);
          
    tTypeIdx = {iBlk.tri.trialType.auditory; iBlk.tri.trialType.visual; iBlk.tri.trialType.coherent; iBlk.tri.trialType.conflict};
    galvoIdx = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3]};
    
    if ~exist('tblData', 'var')
        emCell = cell(length(galvoIdx), length(tTypeIdx));
        tblData = struct('reacT', emCell, 'visC', emCell, 'lasOn', emCell, 'mouseID', emCell);
    end
    
    for j = 1:length(galvoIdx)
        gIdx = ismember(iBlk.tri.inactivation.galvoPosition, galvoIdx{j}, 'rows');
        for k = 1:length(tTypeIdx)
            tBlk = prc.filtBlock(iBlk, (iBlk.tri.inactivation.laserType==0 | gIdx) & tTypeIdx{k});
            normBlk = prc.filtBlock(tBlk, tBlk.tri.inactivation.laserType==0);
            lasBlk = prc.filtBlock(tBlk, tBlk.tri.inactivation.laserType==1);
            lasBlk = prc.filtBlock(lasBlk, lasBlk.tri.stim.visDiff<0 | (lasBlk.tri.stim.visDiff==0 &  lasBlk.tri.stim.audDiff<0));
            
            normGrds = prc.getGridsFromBlock(normBlk, 1);
            lasGrds = prc.getGridsFromBlock(lasBlk, 1);
            
            rN = normGrds.reactionTimeComb;
            rAllN = median(normBlk.tri.outcome.reactionTime);
            rL = lasGrds.reactionTimeComb;
            rAllL = median(lasBlk.tri.outcome.reactionTime);
            visC = abs(normGrds.visValues);
            
            tblData(j,k).reacT = [tblData(j,k).reacT; [rN(~isnan(rN))*0+rAllN; rL(~isnan(rL))*0+rAllL]*1000];
%             tblData(j,k).reacT = [tblData(j,k).reacT; [rN(~isnan(rN)); rL(~isnan(rL))]*1000];
            tblData(j,k).visC = [tblData(j,k).visC; [visC(~isnan(rN)); visC(~isnan(rL))]];
            tblData(j,k).lasOn = [tblData(j,k).lasOn;[rN(~isnan(rN))*0; rL(~isnan(rL))*0+1]];
            tblData(j,k).mouseID = [tblData(j,k).mouseID;[rN(~isnan(rN))*0; rL(~isnan(rL))*0+1]*0+i]; 
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

typeOrd = {'Aud'; 'Vis'; 'Coh'; 'Con'};
lineOpt.Color = [0.5 0.5 0.5];
lineOpt.lineWidth = 1;
for i = 1:2
    for j = 1:4
        axH = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
        hold on
        lIdx = tblData(i,j).lasOn>0;
        
        nData = tblData(i,j).reacT(~lIdx);
        lData = tblData(i,j).reacT(lIdx);
        nMIdx = tblData(i,j).mouseID(~lIdx);
        lMIdx = tblData(i,j).mouseID(lIdx);
        
        meanRTNormLas{j,1} = cell2mat(arrayfun(@(x) [mean(nData(nMIdx==x)) mean(lData(lMIdx==x))], unique(lMIdx), 'uni', 0));
        
        tbl = table(tblData(i,j).reacT,tblData(i,j).lasOn,nominal(tblData(i,j).visC),tblData(i,j).mouseID ...
            ,'VariableNames',{'ReactTime','LaserOn','VisContrast','Mouse'});
        lme = fitlme(tbl, 'ReactTime~LaserOn+VisContrast+(1|Mouse)');
        pVal(j,1) = lme.Coefficients.pValue(strcmpi(lme.Coefficients.Name, 'LaserOn'));
    end
    for j = 1:4
        subMeanRTNormLas = meanRTNormLas{j};
        xDat = ((j-1)*2)+1;
        offset = mean([meanRTNormLas{3}(:,1) meanRTNormLas{4}(:,1)],2);
        subMeanRTNormLas = subMeanRTNormLas-repmat(offset,1,2);
        arrayfun(@(x) plot(axH, [xDat xDat+1], subMeanRTNormLas(x,:), lineOpt), 1:size(subMeanRTNormLas,1))
        plot(axH, [xDat xDat+1], mean(subMeanRTNormLas), 'k', 'linewidth', 2)
    end
    
    
    if i == 1; xlim([0.5 8.5]); ylim([-25 175]); title('MOs inactivation'); end
    if i == 2; xlim([0.5 8.5]); ylim([-35 25]); title('V1 inactivation'); end

    set(gca,'xtick',[1.5:2:7.5],'xticklabel',typeOrd)
    axis square;
    box off;
    legendCell = arrayfun(@(x,y) [x{1} ': p <' num2str(round(y, 2, 'significant'))], typeOrd,pVal,'uni',0);
    legend(legendCell, 'location', 'none', 'FontSize', 10, 'Position', get(gca, 'Position').*[1, 0.7, 1, 0.3])
    legend('boxoff')
end
%%
export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\3_scatterPlotsLMEModel_Shaded_RemovedOutliersNew2', '-pdf', '-painters');
end
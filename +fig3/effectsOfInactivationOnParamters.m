function effectsOfInactivationOnParamters
%%
uniMice = {'PC027'; 'PC029'; 'DJ008'; 'DJ006'; 'DJ007'};
glmsCont_V1_MOs_A1_S1 = cell(length(uniMice), 5);
for i = 1:length(uniMice)
    s = spatialAnalysis(uniMice{i}, 'uniscan',1,1);
    mOpt = struct;
    mOpt.contOnly = 1;
    mOpt.nRepeats = 1;
    mOpt.useDif = 0;
    glmsCont_V1_MOs_A1_S1{i,1} = s.getModelFitsToInactivationData(mOpt);
        
    mOpt.contParams = glmsCont_V1_MOs_A1_S1{i,1}{1}.prmFits;
    mOpt.useDif = 1;
    mOpt.useGroups = 1;
    mOpt.contOnly = 0;
    mOpt.freeP = [1,1,1,0,1,1]>0;
    
    mOpt.groupIDs = 'mos';
    [~, glmsCont_V1_MOs_A1_S1{i,2}] = s.getModelFitsToInactivationData(mOpt);
    mOpt.groupIDs = 'v1';
    [~, glmsCont_V1_MOs_A1_S1{i,3}] = s.getModelFitsToInactivationData(mOpt);
    mOpt.groupIDs = 'a1';
    [~, glmsCont_V1_MOs_A1_S1{i,4}] = s.getModelFitsToInactivationData(mOpt);
    mOpt.groupIDs = 's1';
    [~, glmsCont_V1_MOs_A1_S1{i,5}] = s.getModelFitsToInactivationData(mOpt);
end
%%
tDat = cellfun(@(x) x{1}.prmFits, glmsCont_V1_MOs_A1_S1, 'uni', 0);
tDat(:,1) = [];
clear deltaParams_V1_MOs_A1_S1;
for i = 1:4
    deltaParams_V1_MOs_A1_S1{1,i} = cell2mat(tDat(:,i));
    deltaParams_V1_MOs_A1_S1{1,i} (:,4) = [];
end
delta_mIdx_Prm_Site = cat(3,deltaParams_V1_MOs_A1_S1{:});
prmNames = {'B'; 'VI'; 'VC'; 'AI'; 'AC'};

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


for i = 1:4
    plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg); cla
    yDat = squeeze(delta_mIdx_Prm_Site(:,:,i));
    nXPnts = size(yDat,2);
    xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));
    
    set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
    hold on
    for j = 1:nXPnts-1
        cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(:,j:j+1),2), num2cell(yDat(:,j:j+1),2));
    end
    
    plot(xlim, [0 0], '--k')
    sites ={'MOs'; 'VIS'; 'AUD'; 'S1'};
    plot(gca, xDat, yDat,'ok', 'MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',5);
    [~, pVal] = ttest(yDat);
    xlim([0 nXPnts]);
    set(gca, 'XTick', 0.5:1:4.5, 'XTickLabels', prmNames)
    ylim([-4 2.5]);
    yRng = ylim;
    text(2, yRng(2), sites{i});
    for j = find(pVal < 0.05)
        text(xDat(1,j), max(yDat(:,j))+0.5, '*');
    end
end

%%
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\sigTest4InactivationSites', '-pdf', '-painters');

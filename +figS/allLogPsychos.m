function allLogPsychos(behBlks)
%% This function plots the data panels for figure one of the ms
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1); end
plotIdx = 1;
for i = 1:length(behBlks.blks)
    if mod(plotIdx, 6) == 1
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
    end
    
    singleBlk = spatialAnalysis(behBlks.blks(i).exp.subject{1}, 'behavior', 0, 1);
    singleBlk.blks = prc.filtBlock(singleBlk.blks, singleBlk.blks.tri.stim.visContrast ~= 0.06);
    
    singleBlk.blks = prc.filtBlock(singleBlk.blks, singleBlk.blks.tri.stim.visContrast ~= 0.06);
    axesHandle = plt.tightSubplot(nRows,nCols,i-(floor((i-1)/6)*6),axesGap,botTopMarg,lftRgtMarg);
    singleBlk.viewGLMFits('simpLogSplitVSplitA', [],'log', 1)
    axis square;
    plotIdx = plotIdx+1;
end


%%
export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\SupX_allPsycho3', '-pdf', '-painters'); close
export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\SupX_allPsycho2', '-pdf', '-painters'); close
export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\SupX_allPsycho1', '-pdf', '-painters');
end
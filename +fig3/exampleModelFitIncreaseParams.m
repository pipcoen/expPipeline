function exampleModelFitIncreaseParams(uniBlks)
%%


% %%
% figure;
% axHeight = 250;
% axWidth = 250;
% nCols = 3;
% nRows = 2;
% figHeight = nRows*axHeight;
% figWidth = nCols*axWidth;
% 
% axesGap = [50/figHeight 50/figWidth];
% botTopMarg = [40, 40]/figHeight;
% lftRgtMarg = [40, 40]/figWidth;
% set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);
% 
% for i = 1:3
% fracRightTurns = squeeze(mean(cell2mat(cellfun(@(x) permute(x.fracRightTurns,[3,1,2]), gridsCont_V1_MOs{i}, 'uni', 0)),1));
% fracRightLowB = squeeze(mean(cell2mat(cellfun(@(x) permute(x.fracRightTurnsLowBound,[3,1,2]), gridsCont_V1_MOs{i}, 'uni', 0)),1));
% fracRightHighB = squeeze(mean(cell2mat(cellfun(@(x) permute(x.fracRightTurnsHighBound,[3,1,2]), gridsCont_V1_MOs{i}, 'uni', 0)),1));
% 
% plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg); cla;
% set(gca, 'XTick', [-80 0 80]);
% plotData = cat(3, fracRightTurns, fracRightLowB, fracRightHighB);
% plt.rowsOfGrid(gridsCont_V1_MOs{i}{1}.visValues(1,:)*100, plotData, plt.selectRedBlueColors(gridsCont_V1_MOs{i}{1}.audValues(:,1)));
% end
% %%
% for i = 1:3
%     plt.tightSubplot(nRows,nCols,i+3,axesGap,botTopMarg,lftRgtMarg); cla;
% 
%     allGlMs = glmsCont_V1_MOs{i};
%     glm2Plot = glmsCont_V1_MOs{i}{1};    
%     glm2Plot.prmFits = mean(cell2mat(cellfun(@(x) x.prmFits,allGlMs, 'uni', 0)));
%     glm2Plot.prmInit = mean(cell2mat(cellfun(@(x) x.prmInit,allGlMs, 'uni', 0)));
%     
%     fracRightTurns = squeeze(mean(cell2mat(cellfun(@(x) permute(x.fracRightTurns,[3,1,2]), gridsCont_V1_MOs{i}, 'uni', 0)),1));
%     
%     pHatCalculated = glm2Plot.calculatepHat(glm2Plot.prmFits,'eval');
%     [grids.visValues, grids.audValues] = meshgrid(unique(glm2Plot.evalPoints(:,1)),unique(glm2Plot.evalPoints(:,2)));
%     [~, gridIdx] = ismember(glm2Plot.evalPoints, [grids.visValues(:), grids.audValues(:)], 'rows');
%     plotData = grids.visValues;
%     plotData(gridIdx) = pHatCalculated(:,2);
%     plotOpt.lineStyle = '-';
%     plotOpt.Marker = 'none';
%     
%     contrastPower = 1;
%     visValues = (abs(grids.visValues(1,:))).^contrastPower.*sign(grids.visValues(1,:));
%     lineColors = plt.selectRedBlueColors(grids.audValues(:,1));
%     plt.rowsOfGrid(visValues, plotData, lineColors, plotOpt);
%     
%     plotOpt.lineStyle = 'none';
%     plotOpt.Marker = '.';
%     
%     visValues = gridsCont_V1_MOs{i}{1}.visValues(1,:);
%     maxContrast = max(abs(visValues(1,:)));
%     visValues = abs(visValues).^contrastPower.*sign(visValues)./maxContrast;
%     plt.rowsOfGrid(visValues, fracRightTurns, lineColors, plotOpt);
%    
%     xlim([-1 1])
%     midPoint = 0.5;
%     box off;
%     set(gca, 'xTick', (-1):(1/4):1, 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)*100));
%     
%     xL = xlim; hold on; plot(xL,[midPoint midPoint], '--k', 'linewidth', 1.5);
%     yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
%     
%     currAxPos = get(gca, 'position');
%     axInset = axes;
%     set(gca, 'position', [currAxPos(1)+currAxPos(3)-0.06, currAxPos(2:end).*[1.5 0.25 0.25]])
%     
%     if i == 1
%         X = categorical({'b'; 'VI'; 'VC'; 'Y'; 'AI'; 'AC'});
%         h = bar(X,glm2Plot.prmFits);
%     elseif i == 2
%         X = categorical({'b'; 'VC'});
%         h = bar(X,glm2Plot.prmFits([1 3]));
%     elseif i == 3
%         X = categorical({'b'; 'VI'; 'VC'; 'AI'; 'AC'});
%         h = bar(X,glm2Plot.prmFits([1,2,3,5,6]));
%     end
%     box off;
%     set(h, 'EdgeColor', 'none');
% end

%"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
%see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
% export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_exampleInactModelFits', '-pdf', '-painters');
end
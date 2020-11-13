function visAndM2InactivationBarCharts
%% This function plots the data panels for figure one of the ms
s = spatialAnalysis('all', 'uniscan', 1, 1);

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
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 100 figWidth, figHeight]);
%%
inactivationResults = s.viewInactivationResults('timdif', 1500, {'VL';'AL';'CohL'}, 'avmos'); close;
plt.tightSubplot(nRows,nCols,1:2,axesGap,botTopMarg,lftRgtMarg);
yDat = round(cell2mat(cellfun(@(x) x(~isnan(x))', inactivationResults.meanContData, 'uni', 0))*100);
pVals = cell2mat(cellfun(@(x) x(~isnan(x))', inactivationResults.pVals, 'uni', 0));
%%
pValLabels = cell(size(pVals));
pValLabels(pVals<0.01) = deal({'*'});
pValLabels(pVals<0.001) = deal({'**'});
pValLabels(pVals<0.0001) = deal({'***'});
pValLabels(pVals<0.00001) = deal({'****'});
pValLabels(pVals>0.01) = deal({'NS'});

xDat = 1:length(inactivationResults.subsets);
hB = bar(xDat, yDat);
hT=[];              % placeholder for text object handles
for i=1:length(hB)  % iterate over number of bar objects
  hT=[hT text(hB(i).XData+hB(i).XOffset,hB(i).YData,pValLabels(:,i), ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
end
legend ('AVIpsi', 'M2Ipsi', 'AVContra', 'VisContra')  
set(gca, 'xTickLabel', {'Visual Trials'; 'Auditory Trials'; 'Coherent Trials'})
box off;
ylabel('Change in percentage of timeouts');

%%
inactivationResults = s.viewInactivationResults('readif', 150000, {'CohL', 'ConL'}, 'v1mos'); close;
plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);
yDat = cell2mat(cellfun(@(x) x(1,~isnan(x(1,:))), inactivationResults.meanContData, 'uni', 0))*1000;
pVals = cell2mat(cellfun(@(x) x(1,~isnan(x(1,:))), inactivationResults.pVals, 'uni', 0));

pValLabels = cell(size(pVals));
pValLabels(pVals<0.01) = deal({'*'});
pValLabels(pVals<0.001) = deal({'**'});
pValLabels(pVals<0.0001) = deal({'***'});
pValLabels(pVals<0.00001) = deal({'****'});
pValLabels(pVals>0.01) = deal({'NS'});

xDat = 1:2;
hB = bar(xDat, yDat);
hT=[];              % placeholder for text object handles
for i=1:length(hB)  % iterate over number of bar objects
  hT=[hT text(hB(i).XData+hB(i).XOffset,hB(i).YData,pValLabels(:,i), ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
end

set(gca, 'xTickLabel', {'Coherent Trials'; 'Conflict Trials'})
box off;
ylabel('Change in reaction time (ms)');
legend ('VisIpsi', 'VisContra')  

%%
inactivationResults = s.viewInactivationResults('readif', 15000, {'VL', 'AL', 'CohL', 'ConL'}, 'avmos'); 
close;
plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);
yDat = cell2mat(cellfun(@(x) x(1,~isnan(x(1,:))), inactivationResults.meanContData, 'uni', 0))*1000;
pVals = cell2mat(cellfun(@(x) x(1,~isnan(x(1,:))), inactivationResults.pVals, 'uni', 0));

pValLabels = cell(size(pVals));
pValLabels(pVals<0.01) = deal({'*'});
pValLabels(pVals<0.001) = deal({'**'});
pValLabels(pVals<0.0001) = deal({'***'});
pValLabels(pVals<0.00001) = deal({'****'});
pValLabels(pVals>0.01) = deal({'NS'});

xDat = 1:2;
hB = bar(xDat, yDat);
hT=[];              % placeholder for text object handles
for i=1:length(hB)  % iterate over number of bar objects
  hT=[hT text(hB(i).XData+hB(i).XOffset,hB(i).YData,pValLabels(:,i), ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
end

set(gca, 'xTickLabel', {'Coherent Trials'; 'Conflict Trials'})
box off;
ylabel('Change in reaction time (ms)');
legend ('VisIpsi', 'VisContra')  

%%


export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_modelInactivationBarCharts_AV_100000', '-pdf', '-painters');
end
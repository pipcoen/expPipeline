function inactResultsForRT
%%
load('fig3fInactResultsForRT', 'inactResultsForChoice')
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

subsets = inactResultsForChoice.subsets;
nMice = size(inactResultsForChoice.laserOnData,2);

for i = 1:length(subsets)
    for j = 1:nMice
    nShuffles = size(inactResultsForChoice.shuffLaserOffData{i,j},3);
    contData{i,j} = inactResultsForChoice.meanContEffects{i,j};
    shuffleData = double(inactResultsForChoice.shuffLaserOnData{i,j} - inactResultsForChoice.shuffLaserOffData{i,j});
    sortedData = cellfun(@squeeze, num2cell(sort(abs(cat(3,shuffleData, contData{i,j})),3,'descend'),3), 'uni', 0);
    pVals{i,j} = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(contData{i,j}), sortedData,'uni', 0));
    end
end

contDataV1_MOs = cell2mat(cellfun(@(x) [x(1,5) x(3,4)], contData, 'uni', 0));
pValsV1_MOs = cell2mat(cellfun(@(x) [x(1,5) x(3,4)], pVals, 'uni', 0));

contDataV1_MOs(2,:) = [];
pValsV1_MOs(2,:) = [];

pValLabels = cell(size(pValsV1_MOs));
pValLabels(pValsV1_MOs<0.05) = deal({'*'});
pValLabels(pValsV1_MOs<0.001) = deal({'**'});
pValLabels(pValsV1_MOs<0.0001) = deal({'***'});
pValLabels(pValsV1_MOs<0.00001) = deal({'****'});
pValLabels(pValsV1_MOs>0.05) = deal({'NS'});

xDat = 1:3;
ylim([-25 65]);
plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
hB = bar(xDat, contDataV1_MOs);
hT=[];              % placeholder for text object handles

for i=1:length(hB)  % iterate over number of bar objects
  hT=[hT text(hB(i).XData+hB(i).XOffset,hB(i).YData,pValLabels(:,i), ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
end
legend ('VisContra', 'M2Contra');
set(gca, 'xTickLabel', {'Visual'; 'Coherent'; 'Conflict'});%%
box off;

%"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
%see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_inactResultsForRT', '-pdf', '-painters');
end
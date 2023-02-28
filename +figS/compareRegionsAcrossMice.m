function compareRegionsAcrossMice
% load('fig3aInactResultsForChoice', 'inactResultsForChoice')
load('figS3InactResultsForChoice_IndiMice', 'inactResultsForChoice')
figure;
subsets = inactResultsForChoice.subsets;
nMice = size(inactResultsForChoice.laserOnData,2);
axesOpt.totalNumOfAxes = length(subsets)*nMice;
axesOpt.btlrMargins =  [10 30 10 10];
axesOpt.gapBetweenAxes = [10 0];
axesOpt.axesSize = [200 200];
axesOpt.numOfRows = nMice;
axesOpt.numOfCols = length(subsets);
axesOpt.totalNumOfAxes = length(subsets)*nMice;

subRegions = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4]};
compareEffect_RegSubset = cell(length(subRegions),length(subsets));
%%
gridXY = inactResultsForChoice.gridXY{1};
for j = 1:nMice
    for i = 1:length(subsets)
        for k = 1:length(subRegions)
            contData = inactResultsForChoice.meanContEffects{i,j};
            dat2Take = gridXY{1}*0;
            for q = 1:length(subRegions{k})
                dat2Take(gridXY{1}==subRegions{k}(q,1) & gridXY{2}==subRegions{k}(q,2)) = 1;
            end
            compareEffect_RegSubset{k,i}(j,1) = mean(contData(dat2Take>0));
        end
    end
end

%"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
%see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
% export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\SupX_inactResultsForChoice_IndiMice', '-pdf', '-painters');
end
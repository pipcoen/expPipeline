function viewInactivationResults(obj, plotType)
if ~exist('plotType', 'var'); plotType = 'uni'; end
if ~isempty(obj.expDate)
    switch lower(plotType)
        case {'uni'; 'bil'; 'unisig'}
            runMouseReplicate(copy(obj), {'VL', 'VR', 'AL', 'AR', 'CohL', 'CohR','ConL', 'ConR'}, ['viewInactivationResults(''' lower(plotType) ''')']);
        case {'dif'}
            runMouseReplicate(copy(obj), {'VisUni(L-R)', 'AudUni(L-R)', 'CohUni(VL-VR)', 'ConUni(AL-AR)'}, 'viewInactivationResults(''dif'')');
    end
    return;
end
figure;
axesOpt.totalNumOfAxes = length(obj.subjects);
axesOpt.btlrMargins =  [50 100 10 10];
axesOpt.gapBetweenAxes = [40 0];
axesOpt.figureHWRatio = 0.8;
axesOpt.figureSize = 400;
respBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{1}, 2);
respBlock = prc.combineBlocks(respBlock, respBlock.galvoPosition(:,2)~=4.5);
normBlock = prc.combineBlocks(respBlock, respBlock.laserType==0);
uniBlock = prc.combineBlocks(respBlock, respBlock.laserType==1);
bilBlock = prc.combineBlocks(respBlock, respBlock.laserType==2);

for i  = 1:length(obj.subjects)
    switch lower(plotType)
        case {'uni'; 'bil'}
            if strcmpi(plotType, 'uni'); tempBlock = prc.getDefinedSubset(uniBlock, obj.subjects{i});
            elseif strcmpi(plotType, 'bil'); tempBlock = prc.getDefinedSubset(bilBlock, obj.subjects{i});
                normBlock.galvoPosition(:,1) = abs(normBlock.galvoPosition(:,1));
            end
            
            [scanPlot.data, scanPlot.gridXY] = prc.makeGrid(tempBlock, tempBlock.responseMade==2, @mean, 'galvouni');
            [scanPlot.nTrials] = prc.makeGrid(tempBlock, tempBlock.responseMade==2, @length, 'galvouni');
            
            tempBlock = prc.getDefinedSubset(normBlock, obj.subjects{i});
            scanPlot.data = scanPlot.data - prc.makeGrid(tempBlock, tempBlock.responseMade==2, @mean, 'galvouni');
            axesOpt.numOfRows = 4;
            scanPlot.addTrialNumber = 1;
        case {'unisig'}
            tempBlock = prc.getDefinedSubset(respBlock, obj.subjects{i});
            tempBlock = prc.combineBlocks(tempBlock, tempBlock.timeOutsBeforeResponse==0);
            nShuffles = 10000;
            inactiveGrid = cell(nShuffles+1, 1);
            for j = 1:length(inactiveGrid)
                if j > 1
                    tempBlock.laserType = tempBlock.laserType(randperm(length(tempBlock.laserType)));
                    tempBlock.galvoPosition = tempBlock.galvoPosition(randperm(length(tempBlock.laserType)),:);
                end
                normBlock = prc.combineBlocks(tempBlock, tempBlock.laserType==0);
                laserBlock = prc.combineBlocks(tempBlock, tempBlock.laserType==1);
                [inactiveGrid{j}, scanPlot.gridXY] = prc.makeGrid(laserBlock, laserBlock.responseMade==2, @mean, 'galvouni');
                inactiveGrid{j} = inactiveGrid{j} - mean(normBlock.responseMade==2);
            end
            
            scanPlot.data = reshape(cell2mat(inactiveGrid'), [size(scanPlot.gridXY{1}), nShuffles+1]);
            scanPlot.pVals = nan*scanPlot.data(:,:,1);
            for j = 1:size(scanPlot.data,1)
                for k = 1:size(scanPlot.data,2)
                    if isnan(scanPlot.data(j,k)); continue; end
                    scanPlot.pVals(j,k) = (find(sort(abs(scanPlot.data(j,k,:)), 'descend')==abs(scanPlot.data(j,k,1)),1)-1)./nShuffles;
                end
            end
            
            scanPlot.data = scanPlot.data(:,:,1);
            scanPlot.addTrialNumber = 0;
            axesOpt.numOfRows = 4;
        case {'dif'}
            runMouseReplicate(copy(obj), {'VisUni(L-R)', 'AudUni(L-R)', 'CohUni(VL-VR)', 'ConUni(AL-AR)'}, 'viewInactivationResults(''dif'')');
            [scanPlot.plotData, scanPlot.pVals, scanPlot.gridXY, scanPlot.nTrials] = fit.testLaserEffectAtEachSite(respBlock, obj.subjects{i});
    end
    obj.axesHandles = plt.getAxes(axesOpt, i);
    scanPlot.title = obj.subjects{i};
    plt.scanningBrainEffects(scanPlot);
end
end
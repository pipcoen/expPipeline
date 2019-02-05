function viewInactivationResults(obj, plotType)
if ~exist('plotType', 'var'); plotType = 'uni'; end
if ~isempty(obj.expDate)
    switch lower(plotType)
        case {'uni'; 'bil'; 'unisig'; 'unisigrea'}
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
respBlock = prc.filtStruct(obj.blocks{1}, obj.blocks{1}.responseMade~=0);
respBlock = prc.filtStruct(respBlock, respBlock.galvoPosition(:,2)~=4.5);
respBlock = prc.filtStruct(respBlock, ~ismember(abs(respBlock.galvoPosition(:,1)),[0.5; 2; 3.5; 5]));
normBlock = prc.filtStruct(respBlock, respBlock.laserType==0);
uniBlock = prc.filtStruct(respBlock, respBlock.laserType==1);
bilBlock = prc.filtStruct(respBlock, respBlock.laserType==2);
if ~isempty(bilBlock); bilBlock.galvoPosition(:,1) = abs(bilBlock.galvoPosition(:,1)); end

for i  = 1:length(obj.subjects)
    switch lower(plotType)
        case {'uni'; 'bil'}
            if strcmpi(plotType, 'uni')
                tempBlock = prc.getDefinedSubset(uniBlock, obj.subjects{i});
            elseif strcmpi(plotType, 'bil')
                tempBlock = prc.getDefinedSubset(bilBlock, obj.subjects{i});
                normBlock.galvoPosition(:,1) = abs(normBlock.galvoPosition(:,1));
            end
            tempBlock = prc.filtStruct(tempBlock, tempBlock.timeOutsBeforeResponse==0);
            [nonUniformData, scanPlot.gridXY] = prc.makeGrid(tempBlock, [tempBlock.responseMade==2 tempBlock.subjectIdx], [], 'galvouni',2);
            scanPlot.data = cellfun(@(x) mean(x(prc.makeFreqUniform(x(:,2)),1)), nonUniformData);
            scanPlot.nTrials = cellfun(@(x) length(x(prc.makeFreqUniform(x(:,2)),1)), nonUniformData);
 
            tempBlock = prc.getDefinedSubset(normBlock, obj.subjects{i});
            controlData = prc.makeGrid(tempBlock, [tempBlock.responseMade==2 tempBlock.subjectIdx], [], 'galvouni',2);
            controlData = cellfun(@(x) mean(x(prc.makeFreqUniform(x(:,2)),1)), controlData);
            
            scanPlot.data = scanPlot.data - controlData;
            axesOpt.numOfRows = 4;
            scanPlot.addTrialNumber = 1;
        case {'unisig'; 'unisigrea'}
            tempBlock = prc.getDefinedSubset(respBlock, obj.subjects{i});
            tempBlock = prc.filtStruct(tempBlock, tempBlock.timeOutsBeforeResponse==0);  
            nShuffles = 100;
            normEstRepeats = round(nShuffles/10);
            inactiveGrid = cell(nShuffles+1, 1);
            for j = 1:(nShuffles + normEstRepeats)
                laserBlock = prc.filtStruct(tempBlock, tempBlock.laserType==1);
                [nonUniformLaser] = prc.makeGrid(laserBlock, [(1:length(laserBlock.subjectIdx))' laserBlock.subjectIdx], [], 'galvouni',2);
                uniformLaser = cellfun(@(x) x(prc.makeFreqUniform(x(:,2)),1), nonUniformLaser, 'uni', 0);
                uniformLaserLogical = ismember((1:length(laserBlock.subjectIdx))',cell2mat(uniformLaser(:))); 
                laserBlock = prc.filtStruct(laserBlock, uniformLaserLogical>0);
                
                normBlock = prc.filtStruct(tempBlock, tempBlock.laserType==0);
                normBlock = prc.filtStruct(normBlock, prc.makeFreqUniform(normBlock.subjectIdx));
                
                subsampledBlock = prc.combineBlocks([laserBlock;normBlock]);
                if j > normEstRepeats
                    subsampledBlock.laserType = tempBlock.laserType(randperm(length(subsampledBlock.laserType)));
                    subsampledBlock.galvoPosition = tempBlock.galvoPosition(randperm(length(subsampledBlock.laserType)),:);
                end
                
                normBlock = prc.filtStruct(subsampledBlock, subsampledBlock.laserType==0);
                laserBlock = prc.filtStruct(subsampledBlock, subsampledBlock.laserType==1);
                if strcmp(lower(plotType), 'unisigrea')
                    [inactiveGrid{j}, scanPlot.gridXY] = prc.makeGrid(laserBlock, laserBlock.timeToWheelMove, @median, 'galvouni');
                else
                    [inactiveGrid{j}, scanPlot.gridXY] = prc.makeGrid(laserBlock, laserBlock.responseMade==2, @mean, 'galvouni');
                end
                inactiveGrid{j} = inactiveGrid{j} - mean(normBlock.responseMade==2);
            end
            %%
            trueData = mean(reshape(cell2mat(inactiveGrid(1:normEstRepeats)'), [size(scanPlot.gridXY{1}), normEstRepeats]),3);
            shuffleData = reshape(cell2mat(inactiveGrid(normEstRepeats+1:end)'), [size(scanPlot.gridXY{1}), nShuffles]);
            shuffleData = cellfun(@squeeze, num2cell(sort(cat(3,shuffleData, trueData),3,'descend'),3), 'uni', 0);
            scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, trueData, shuffleData,'uni', 0));

            scanPlot.data = trueData;
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
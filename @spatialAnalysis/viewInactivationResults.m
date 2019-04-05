function viewInactivationResults(obj, plotType, nShuffles)
if ~exist('plotType', 'var'); plotType = 'unires'; end
if ~exist('nShuffles', 'var'); nShuffles = 0; end
plotType = lower(plotType);
if ~isempty(obj.expDate)
    if contains(plotType, 'dif')
        subsets = {'VL', 'AL', 'CohL','ConL'};
    else
        subsets = {'VL', 'VR', 'AL', 'AR', 'CohL', 'CohR','ConL', 'ConR'};
    end
    runMouseReplicate(copy(obj), subsets, ['viewInactivationResults(''' plotType ''',' num2str(nShuffles) ')']);
    return;
end
figure;
axesOpt.totalNumOfAxes = length(obj.subjects);
axesOpt.btlrMargins =  [50 100 10 10];
axesOpt.gapBetweenAxes = [40 0];
axesOpt.figureHWRatio = 0.8;
axesOpt.figureSize = 400;
%%
initBlock = prc.filtStruct(obj.blocks{1}, obj.blocks{1}.galvoPosition(:,2)~=4.5);
initBlock = prc.filtStruct(initBlock, ~ismember(abs(initBlock.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | initBlock.laserType==0);
initBlock = prc.filtStruct(initBlock, initBlock.timeOutsBeforeResponse==0);

op2Use = @mean;

subjectIndexes = unique(initBlock.subjectIdx);
condLabels = unique(initBlock.conditionLabelRow(:,1));

if contains(plotType, 'grp')
    galvoGrps = {[1.8, -3; 1.8 -4; 1.8 -2; 3, -4; 3, -3; 3, -2; 4.2, -4; 4.2, -3; 4.2, -2;];[0.6 2; 1.8, 2; 0.6, 3]};%;[4.2 -2; 4.2, -3]};
%     galvoGrps = {[1.8 -4; 3, -4; 3, -3];[0.6 2; 1.8, 2; 0.6, 3]};%;[4.2 -2; 4.2, -3]};
    galvoGrps = [galvoGrps; cellfun(@(x) [-x(:,1) x(:,2)], galvoGrps, 'uni', 0)];
    groupedBlock = cellfun(@(x) prc.filtStruct(initBlock, ismember(initBlock.galvoPosition, x, 'rows') & initBlock.laserType~=0), galvoGrps);
    meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
    meanPositions = cellfun(@(x,y) x*0+repmat(y, size(x,1),1), {groupedBlock.galvoPosition}', meanPositions, 'uni', 0);
    [groupedBlock.galvoPosition] = deal(meanPositions{:});
    initBlock = prc.combineBlocks([groupedBlock; prc.filtStruct(initBlock, initBlock.laserType==0)]);
end

if contains(plotType, 'dif')
    initBlock = prc.switchBlockToContraIpsi(initBlock);
end
%%
if contains(plotType, 'res')
    initBlock = prc.filtStruct(initBlock, initBlock.responseMade~=0);
    initBlock.data2Use = initBlock.responseMade==2;
    scanPlot.colorBarLimits = [-0.8 0.8];
elseif contains(plotType, 'tim')
    initBlock.data2Use = (initBlock.responseMade==0);
    scanPlot.colorBarLimits = [-0.25 0.25];
elseif contains(plotType, 'rea')
    initBlock = prc.filtStruct(initBlock, initBlock.responseMade~=0);
    reactionTimes = initBlock.timeToWheelMove;
    reactionTimes(initBlock.laserType~=0) = nan;
    [tDat1,tDat2] = meshgrid(subjectIndexes,condLabels);
    allCombos = num2cell([tDat2(:) tDat1(:)],2);
    refData = [initBlock.conditionLabelRow(:,1) initBlock.subjectIdx];
    normReactionTimes = cellfun(@(x) nanmedian(reactionTimes(all(refData==x,2))),allCombos);
    [~,normReactionIdx] = ismember(refData, cell2mat(allCombos), 'rows');
    initBlock.data2Use = initBlock.timeToWheelMove./normReactionTimes(normReactionIdx);
    op2Use = @median;
    scanPlot.colorBarLimits = [-0.2 0.2];
end
normBlock = prc.filtStruct(initBlock, initBlock.laserType==0);
uniBlock = prc.filtStruct(initBlock, initBlock.laserType==1);
bilBlock = prc.filtStruct(initBlock, initBlock.laserType==2);

for i  = 1:length(obj.subjects)
    if ~contains(plotType, {'sig'})
        if ~nShuffles; nShuffles = 50; end
        contBlock = prc.getDefinedSubset(normBlock, obj.subjects{i});
        contData = repmat(contBlock.data2Use,1,nShuffles);
        subSamples = prc.makeFreqUniform(contBlock.subjectIdx,nShuffles);
        contData(~subSamples) = nan;
        contData = mean(nanmedian(contData));
        
        if strcmp(plotType(1:3), 'uni'); testBlock = prc.getDefinedSubset(uniBlock, obj.subjects{i});
        elseif strcmp(plotType(1:3), 'bil'); testBlock = prc.getDefinedSubset(bilBlock, obj.subjects{i});
        end
        
        trialIdx = (1:length(testBlock.subjectIdx))';
        nonUniformLaser = prc.makeGrid(testBlock, [testBlock.subjectIdx trialIdx], [], 'galvouni',2);
        laserShuffles = cellfun(@(x) prc.makeFreqUniform(x(:,1),nShuffles,x(:,2)), nonUniformLaser, 'uni', 0);
        laserShuffles = num2cell(cell2mat(laserShuffles(:)),1);
        uniformLaserFilters = cellfun(@(x) ismember(trialIdx, x), laserShuffles, 'uni', 0);
        
        inactiveGrid = cell(nShuffles,1);
        for j =1:nShuffles
            subTestBlock = prc.filtStruct(testBlock, uniformLaserFilters{j});
            [inactiveGrid{j}, scanPlot.gridXY] = prc.makeGrid(subTestBlock, subTestBlock.data2Use, op2Use, 'galvouni');
        end
        scanPlot.data = mean(cat(3,inactiveGrid{:}),3)-contData;
        scanPlot.nTrials = prc.makeGrid(subTestBlock, subTestBlock.data2Use, @length, 'galvouni');
        
        axesOpt.numOfRows = 4;
        scanPlot.addTrialNumber = 1;
    elseif contains(plotType, {'sig'})
        if ~nShuffles; nShuffles = 15000; end
        normEstRepeats = round(nShuffles/10);
        inactiveGrid = cell(nShuffles+normEstRepeats, 1);
        totalLoops = nShuffles + normEstRepeats;
        testBlock = prc.getDefinedSubset(uniBlock, obj.subjects{i});
        contBlock = prc.getDefinedSubset(normBlock, obj.subjects{i});
        
        contBlock.galvoPosition = testBlock.galvoPosition(randi(size(testBlock.galvoPosition,1), size(contBlock.galvoPosition,1),1),:);
        %%
        trialIdx = (1:length(testBlock.subjectIdx))';
        nonUniformLaser = prc.makeGrid(testBlock, [testBlock.subjectIdx trialIdx], [], 'galvouni',2);
        laserShuffles = cellfun(@(x) prc.makeFreqUniform(x(:,1),totalLoops,x(:,2)), nonUniformLaser, 'uni', 0);
        laserShuffles = num2cell(cell2mat(laserShuffles(:)),1);
        uniformLaserFilters = cellfun(@(x) ismember(trialIdx, x), laserShuffles, 'uni', 0);
        
        for j = 1:(totalLoops)
            subTestBlock = prc.filtStruct(testBlock, uniformLaserFilters{j});
            subContBlock = prc.filtStruct(contBlock, prc.makeFreqUniform(contBlock.subjectIdx));
            
            randomBlock = prc.combineBlocks([subTestBlock;subContBlock]);
            if j > normEstRepeats
                randomBlock.laserType = randomBlock.laserType(randperm(length(randomBlock.laserType)));
                randomBlock.galvoPosition = randomBlock.galvoPosition(randperm(length(randomBlock.laserType)),:);
            end
            
            subContBlock = prc.filtStruct(randomBlock, randomBlock.laserType==0);
            subTestBlock = prc.filtStruct(randomBlock, randomBlock.laserType==1);
            [inactiveGrid{j}, scanPlot.gridXY] = prc.makeGrid(subTestBlock, subTestBlock.data2Use, op2Use, 'galvouni');
            inactiveGrid{j} = inactiveGrid{j} - op2Use(subContBlock.data2Use);
        end
        %%
        contData = mean(cat(3,inactiveGrid{1:normEstRepeats}),3);
        shuffleData = cat(3,inactiveGrid{normEstRepeats+1:end});
        sortedData = cellfun(@squeeze, num2cell(sort(abs(cat(3,shuffleData, contData)),3,'descend'),3), 'uni', 0);
        scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(contData), sortedData,'uni', 0));
        
        scanPlot.data = contData; disp(contData);
        scanPlot.addTrialNumber = 0;
        axesOpt.numOfRows = 4;
    end
    obj.hand.axes = plt.getAxes(axesOpt, i);
    scanPlot.title = obj.subjects{i};
    plt.scanningBrainEffects(scanPlot);
end
end
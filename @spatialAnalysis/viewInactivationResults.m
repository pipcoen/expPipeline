function viewInactivationResults(obj, plotType, nShuffles, subsets)
%% Method for "spatialAnalysis" class. Plots effects of inactivation on behavior. Plots are shown as grids on an outline of cortex.

%INPUTS(default values)
%plotType(res)---------A string that contains three letter tags indicates the type of plot to generate. Can be combination of below options
%   'res'-------------------------quantify changes in the mouse response (i.e. fraciton of rightward choices)
%   'dif'-------------------------rather than separately analysing left and right trials, combine trials and use ipsilateral and contralateral
%   'grp'-------------------------combine inactivation sites into M2 and Vis
%   'sig'-------------------------test the significance of inactivation by shuffling inactivation sites and laser on/off
%nShuffles-------------The number of times to shuffle te data
%subsets---------------The data subsets to make plots for (see function "prc.getDefinedSubset")

%Set up defaults for the input values. "op2use" is "mean" as default for responses, but changes below depending on the type data being used
plotType = lower(plotType);
if ~exist('plotType', 'var'); plotType = 'res'; end
if ~exist('nShuffles', 'var'); nShuffles = 0; end
if ~contains(plotType, 'dif') && ~exist('subsets', 'var');subsets = {'VL', 'VR', 'AL', 'AR', 'CohL', 'CohR','ConL', 'ConR'}; end
if contains(plotType, 'dif') && ~exist('subsets', 'var'); subsets = {'VL', 'AL', 'CohL', 'ConL'}; end
if ~iscell(subsets); subsets = {subsets}; end
numOfMice = length(obj.blks);

%Set up plotting axes arrangement on figrure
figure;
axesOpt.totalNumOfAxes = length(subsets)*numOfMice;
axesOpt.btlrMargins =  [10 30 10 10];
axesOpt.gapBetweenAxes = [10 0];
axesOpt.axesSize = [200 200];

for mIdx = 1:numOfMice
    %Create "iBlk" (initialBlock) which removes some incorrect galvoPositions, repeated trials, and keeps only valid trials
    iBlk = prc.filtBlock(obj.blks(mIdx), obj.blks(mIdx).tri.inactivation.galvoPosition(:,2)~=4.5 | obj.blks(mIdx).tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
    
    %Conditional to optionally group inactivaiton sites together (e.g. if you want to combine all V1 sites). We "reflect" these groups, so only one
    %hemisphere needs to be defined. We find all trials with galvoPositions belonging to those groups in iBlk and replace with the mean group position
    if contains(plotType, 'grp')
        %         galvoGrps = {[1.8 -4; 3,-4; 3,-3; 1.8,-3; 4.2,-4; 4.2,-3; 1.8,-2; 3,-2; 4.2,-2];[0.6 2; 1.8, 2; 0.6, 3]};
                galvoGrps = {[1.8 -4; 3,-4; 3,-3;]; [4.2,-3; 4.2,-2];[0.6 2; 1.8, 2; 0.6, 3]};
%         galvoGrps = {[1.8 -4; 3,-4; 3,-3];[0.6 2; 1.8, 2; 0.6, 3]};
        galvoGrps = [galvoGrps; cellfun(@(x) [x(:,1)*-1, x(:,2)], galvoGrps, 'uni', 0)];
        grpNum = length(galvoGrps);
        grpIdx = cellfun(@(x,y) ismember(iBlk.tri.inactivation.galvoPosition, x, 'rows').*y, galvoGrps, num2cell(1:grpNum)', 'uni', 0);
        grpIdx = sum(cell2mat(grpIdx'),2);
        meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
        iBlk.tri.inactivation.galvoPosition(grpIdx>0,:) = cell2mat(meanPositions(grpIdx(grpIdx>0)));
        iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | grpIdx>0);
    end
    
    %If plotType contains 'dif' then we want to switch the responseMade, vis and aud paramters such that all "left" trials are flipped (vis left trials in
    %the case of conflict trials). Now, inactivations on the right hemisphere are contralateral, and left hemisphere is ipsilateral
    if contains(plotType, 'dif')
        rIdx = iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0);
        iBlk.tri.outcome.responseMade(rIdx) = (iBlk.tri.outcome.responseMade(rIdx)*-1+3).*(iBlk.tri.outcome.responseMade(rIdx)>0);
        iBlk.tri.stim.audInitialAzimuth(rIdx) = iBlk.tri.stim.audInitialAzimuth(rIdx)*-1;
        iBlk.tri.stim.visInitialAzimuth(rIdx) = iBlk.tri.stim.visInitialAzimuth(rIdx)*-1;
        iBlk.tri.stim.visInitialAzimuth(isinf(iBlk.tri.stim.visInitialAzimuth)) = inf;
        iBlk.tri.inactivation.galvoPosition(rIdx,1) = -1*iBlk.tri.inactivation.galvoPosition(rIdx,1);
    end
    
    %We add "data2Use" to the outcome structure, and this will change depend on the data we want to compare between inactivaiton sites
    if contains(plotType, 'res')
        %Changes in the fraction of rightward choices
        iBlk = prc.filtBlock(iBlk, iBlk.tri.outcome.responseMade~=0);
        iBlk.tri.inactivation.data2Use = iBlk.tri.outcome.responseMade==2;
        scanPlot.colorBarLimits = [-0.7 0.7];
        op2Use = @mean;
    elseif contains(plotType, 'tim')
        %Changes in the number of timeout trials
        iBlk.tri.inactivation.data2Use = (iBlk.tri.outcome.responseMade==0);
        scanPlot.colorBarLimits = [-0.25 0.25];
        op2Use = @mean;
    elseif contains(plotType, 'rea')
        %Changes in the median reaction time
        iBlk = prc.filtBlock(iBlk, iBlk.tri.outcome.responseMade~=0 & ~isnan(iBlk.tri.outcome.threshMoveTime));
        data2Use = iBlk.tri.outcome.threshMoveTime;
        scanPlot.colorBarLimits = [-0.01 0.01];
        for i = 1:iBlk.tot.subjects
            subIdx = iBlk.tri.subjectRef==i;
%             data2Use(subIdx) = (data2Use(subIdx)-median(data2Use(subIdx)))./mad(data2Use(subIdx),1);
            data2Use(subIdx) = (data2Use(subIdx)/median(data2Use(subIdx)));
            scanPlot.colorBarLimits = [-0.1 0.1];
        end
        iBlk.tri.inactivation.data2Use = data2Use;
        op2Use = @median;
    end
    
    %Create normBlk and uniBlk which are filtered versions of iBlk with only control or inactivation trials respectively
    normBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0);
    uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==1);
    
    for i  = 1:length(subsets)
        %Get subset of trials based on tag. e.g. 'vl' for visual left for both normBlk (contBlk) and uniBlk (testBlk). We need these irrespective of 'sig'
        contBlk = prc.getDefinedSubset(normBlk, subsets{i});
        testBlk = prc.getDefinedSubset(uniBlk, subsets{i});
               
        %Define number of suffles, and the number of times to estimate the control (default is 10% or 500) since this will change with different
        %subsamples of the subjetcs. Total loops is the sum of shuffle and control loops. We randomly assign galvoPositions to the contBlk from the
        %galvoPositions in the testBlock (for shuffling purposes)
        if ~nShuffles; nShuffles = 15000; end
        normEstRepeats = min([round(nShuffles/10) 500]);
        totalLoops = nShuffles + normEstRepeats;
        contBlk.tri.inactivation.galvoPosition = testBlk.tri.inactivation.galvoPosition(randi(testBlk.tot.trials, contBlk.tot.trials,1),:);
        
        %Use some matrix tricks to create "uniformLaserFilters" which is a filter for each shuffle that equalizes the frequency of subjects
        %contributing to each point in the grid of galvo positions.
        trialIdx = (1:testBlk.tot.trials)';
        nonUniformLaser = prc.makeGrid(testBlk, [double(testBlk.tri.subjectRef) trialIdx], [], 'galvouni',2);
        laserShuffles = cellfun(@(x) double(prc.makeFreqUniform(x(:,1),totalLoops,x(:,2))), nonUniformLaser, 'uni', 0);
        laserShuffles = num2cell(cell2mat(laserShuffles(:)),1);
        uniformLaserFilters = cellfun(@(x) ismember(trialIdx, x), laserShuffles, 'uni', 0);
        
        %Removing these excess fields makes "filtBlock" run significantly faster in the subsequent loop
        testBlk.tri = rmfield(testBlk.tri, {'trialType', 'timings', 'outcome', 'stim'});
        contBlk.tri = rmfield(contBlk.tri, {'trialType', 'timings', 'outcome', 'stim'});
        testBlk.tri.inactivation = rmfield(testBlk.tri.inactivation, {'laserPower', 'galvoType', 'laserOnsetDelay', 'laserDuration'});
        contBlk.tri.inactivation = rmfield(contBlk.tri.inactivation, {'laserPower', 'galvoType', 'laserOnsetDelay', 'laserDuration'});
        
        %This is the same loop as above, to generate "inactiveGrid" but with some extra steps to deal with the shuffles
        inactiveGrid = cell(nShuffles+normEstRepeats, 1);
        for j = 1:(totalLoops)
            %Filter both blks to make contributions from each subject equal at each galvoPosition and across control trials.
            subTestBlk = prc.filtBlock(testBlk, uniformLaserFilters{j});
            subContBlk = prc.filtBlock(contBlk, prc.makeFreqUniform(contBlk.tri.subjectRef));
            subContBlk.tri.inactivation.galvoPosition = subTestBlk.tri.inactivation.galvoPosition(randi(subTestBlk.tot.trials, subContBlk.tot.trials,1),:);
            
            %Generate "randomBlk" by concatenating control and test blocks. If the repeat number is outside the "normEstRepeats" range then we shuffle
            %the laser activation and galvoPositions
            randomBlk = prc.catStructs([subTestBlk;subContBlk]);
            if j > normEstRepeats
                randomBlk.tri.inactivation.laserType = randomBlk.tri.inactivation.laserType(randperm(randomBlk.tot.trials));
                randomBlk.tri.inactivation.galvoPosition = randomBlk.tri.inactivation.galvoPosition(randperm(randomBlk.tot.trials),:);
            end
            
            %After shuffling, separate randomBlk based on laser activation. This will do nothing if unshuffled of course. Then estimate the result
            %(which is based on data2Use and op2Use) for each galvoPosition for subTestBlk and subtract the single value from the subContBlk.
            subContBlk = prc.filtBlock(randomBlk, randomBlk.tri.inactivation.laserType==0);
            subTestBlk = prc.filtBlock(randomBlk, randomBlk.tri.inactivation.laserType==1);
            [inactiveGrid{j}, scanPlot.gridXY] = prc.makeGrid(subTestBlk, subTestBlk.tri.inactivation.data2Use, op2Use, 'galvouni',[],1);
            inactiveGrid{j} = inactiveGrid{j} - op2Use(subContBlk.tri.inactivation.data2Use);
        end
        
        %"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
        %see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
        contData = mean(cat(3,inactiveGrid{1:normEstRepeats}),3);
        shuffleData = cat(3,inactiveGrid{normEstRepeats+1:end});
        sortedData = cellfun(@squeeze, num2cell(sort(abs(cat(3,shuffleData, contData)),3,'descend'),3), 'uni', 0);
        scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(contData), sortedData,'uni', 0));
        scanPlot.data = contData; disp(contData);
        scanPlot.addTrialNumber = 0;
        sigLevels = (10.^(-2:-1:-10))';
        lastSigLevel = find(sigLevels>min(scanPlot.pVals(:)),1,'last');
        scanPlot.sigLevels = sigLevels(max([1 lastSigLevel-2]):lastSigLevel);
        scanPlot.sigLevels = [0.05; 0.01; 0.001];
        
        %%Plot the data in the "scanPlot" structure.
        axesOpt.numOfRows = 1;%min([4 length(subsets)]);
        axesOpt.numOfCols = 4;%numOfMice;
        obj.hand.axes = plt.getAxes(axesOpt, ((i-1)*numOfMice)+mIdx);
        scanPlot.title = subsets{i};
        plt.scanningBrainEffects(scanPlot);
    end
end
end
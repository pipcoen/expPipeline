function inactResultsForModel = viewInactivationEffectsOnModel(obj, plotType, nShuffles, freeP, groups)
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
numOfMice = length(obj.blks);
if numOfMice > 1; error('Only coded to handle one mouse atm'); end
if ~exist('groups', 'var'); useGroups = 0; else; useGroups = 1; end
if ~exist('nShuffles', 'var'); nShuffles = 0; end
if ~exist('plotType', 'var'); plotType = 'grp'; end
if ~exist('freeP', 'var'); freeP = [1 1 1 0 1 1]>0; end
freeP = freeP>0;

%Create "iBlk" (initialBlock) which removes some incorrect galvoPositions, repeated trials, and keeps only valid trials
iBlk = prc.filtBlock(obj.blks, obj.blks.tri.inactivation.galvoPosition(:,2)~=4.5);
iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
iBlk = prc.filtBlock(iBlk, iBlk.tri.stim.visContrast~=0.06);

%Conditional to optionally group inactivaiton sites together (e.g. if you want to combine all V1 sites). We "reflect" these groups, so only one
%hemisphere needs to be defined. We find all trials with galvoPositions belonging to those groups in iBlk and replace with the mean group position
if useGroups
    galvoGrps = {};
    groups = lower(groups);
    if contains(groups, 'v1'); galvoGrps = [galvoGrps; [1.8 -4; 3,-4; 3,-3]]; end
    if contains(groups, 'a1'); galvoGrps = [galvoGrps; [4.2,-2; 4.2,-3; 4.2,-4]]; end
    if contains(groups, 'mos'); galvoGrps = [galvoGrps; [0.6 2; 1.8, 2; 0.6, 3]]; end
    if contains(groups, 'av'); galvoGrps = [galvoGrps; [1.8 -4; 3,-4; 4.2,-4; 1.8,-3; 3,-3; 4.2,-3; 1.8,-2; 3,-2; 4.2,-2]]; end
    
    galvoGrps = [galvoGrps; cellfun(@(x) [x(:,1)*-1, x(:,2)], galvoGrps, 'uni', 0)];
    grpNum = length(galvoGrps);
    grpIdx = cellfun(@(x,y) ismember(iBlk.tri.inactivation.galvoPosition, x, 'rows').*y, galvoGrps, num2cell(1:grpNum)', 'uni', 0);
    grpIdx = sum(cell2mat(grpIdx'),2);
    meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
    iBlk.tri.inactivation.galvoPosition(grpIdx>0,:) = cell2mat(meanPositions(grpIdx(grpIdx>0)));
    iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | grpIdx>0);
end

contBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 0);
uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 1);

if contains(plotType, 'dif')
    %If plotType contains 'dif' then we want to switch the responseCalc, vis and aud paramters such that all "left" trials are flipped (vis left trials in
    %the case of conflict trials). Now, inactivations on the right hemisphere are contralateral, and left hemisphere is ipsilateral
    idx2Flip = uniBlk.tri.inactivation.galvoPosition(:,1)<0;
    uniBlk.tri.outcome.responseCalc(idx2Flip) = (uniBlk.tri.outcome.responseCalc(idx2Flip)*-1+3).*(uniBlk.tri.outcome.responseCalc(idx2Flip)>0);
    uniBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*uniBlk.tri.inactivation.galvoPosition(idx2Flip,1);
    uniBlk.tri.stim.audDiff(idx2Flip) = -1*uniBlk.tri.stim.audDiff(idx2Flip);
    uniBlk.tri.stim.visDiff(idx2Flip) = -1*uniBlk.tri.stim.visDiff(idx2Flip);
    uniBlk.tri.stim.conditionLabel(idx2Flip) = -1*uniBlk.tri.stim.conditionLabel(idx2Flip);
end

[~, gridXY] = prc.makeGrid(uniBlk, uniBlk.tri.outcome.responseCalc, [], 'galvouni',2);

%Define number of suffles, and the number of times to estimate the control (default is 10% or 500) since this will change with different
%subsamples of the subjetcs. Total loops is the sum of shuffle and control loops. We randomly assign galvoPositions to the contBlk from the
%galvoPositions in the testBlock (for shuffling purposes)
if ~nShuffles; nShuffles = 1500; end
nShuffles = round(nShuffles/10)*10;
normEstRepeats = round(nShuffles/10);

%Use some matrix tricks to create "uniformLaserFilters" which is a filter for each shuffle that equalizes the frequency of subjects
%contributing to each point in the grid of galvo positions.
trialIdx = (1:uniBlk.tot.trials)';
nonUniformLaser = prc.makeGrid(uniBlk, [double(uniBlk.tri.subjectRef) trialIdx], [], 'galvouni',2);
laserShuffles = cellfun(@(x) double(prc.makeFreqUniform(x(:,1),normEstRepeats,x(:,2))), nonUniformLaser, 'uni', 0);
laserShuffles = num2cell(cell2mat(laserShuffles(:)),1)';
uniformLaserFilters = repmat(cellfun(@(x) ismember(trialIdx, x), laserShuffles, 'uni', 0),11,1);
uniformControlFilters = repmat(num2cell(prc.makeFreqUniform(contBlk.tri.subjectRef,normEstRepeats),1)',11,1);

exampleUniform = prc.filtBlock(uniBlk, uniformLaserFilters{1});
totalLaserTrials = prc.makeGrid(exampleUniform, exampleUniform.tri.outcome.responseCalc, @length, 'galvouni');

%Removing these excess fields makes "filtBlock" run significantly faster in the subsequent loop
uniBlk.tri = rmfield(uniBlk.tri, {'timings'});
contBlk.tri = rmfield(contBlk.tri, {'timings'});
uniBlk.tri.inactivation = rmfield(uniBlk.tri.inactivation, {'laserPower', 'galvoType', 'laserOnsetDelay', 'laserDuration'});
contBlk.tri.inactivation = rmfield(contBlk.tri.inactivation, {'laserPower', 'galvoType', 'laserOnsetDelay', 'laserDuration'});

%This is the same loop as above, to generate "inactiveGrid" but with some extra steps to deal with the shuffles
[contParams, deltaParams] = deal(cell(size(gridXY{1})));
for i = 1:normEstRepeats+nShuffles
    %Filter both blks to make contributions from each subject equal at each galvoPosition and across control trials.
    uniformLasBlk = prc.filtBlock(uniBlk, uniformLaserFilters{i});
    uniformContBlk = prc.filtBlock(contBlk, uniformControlFilters{i});
    uniformContBlk.tri.inactivation.galvoPosition = uniformLasBlk.tri.inactivation.galvoPosition(randi(uniformLasBlk.tot.trials, uniformContBlk.tot.trials,1),:);

    %Generate "randomBlk" by concatenating control and test blocks. If the repeat number is outside the "normEstRepeats" range then we shuffle
    %the laser activation and galvoPositions
    randomBlk = prc.catStructs([uniformLasBlk;uniformContBlk]);
    if i > normEstRepeats
        randomBlk.tri.inactivation.laserType = randomBlk.tri.inactivation.laserType(randperm(randomBlk.tot.trials));
        randomBlk.tri.inactivation.galvoPosition = randomBlk.tri.inactivation.galvoPosition(randperm(randomBlk.tot.trials),:);
    end
    subContBlk = prc.filtBlock(randomBlk, randomBlk.tri.inactivation.laserType==0);
    subTestBlk = prc.filtBlock(randomBlk, randomBlk.tri.inactivation.laserType==1);
    
    contGLM = fit.GLMmulti(subContBlk, 'simpLogSplitVSplitA');
    contGLM.fit;
    
    for j = 1:length(gridXY{1}(:))
        if totalLaserTrials(j) ==0; deltaParams{j}(i,:) = nan*ones(1,6); continue; end
        contParams{j}(i,:) = contGLM.prmFits;
        galvoRef = [gridXY{1}(j) gridXY{2}(j)];
        includeIdx = ismember(subTestBlk.tri.inactivation.galvoPosition,galvoRef, 'rows');
        subTestGalvoBlk = prc.filtBlock(subTestBlk, includeIdx);
        
        %After shuffling, separate randomBlk based on laser activation. This will do nothing if unshuffled of course. Then estimate the result
        %(which is based on data2Use and op2Use) for each galvoPosition for subTestBlk and subtract the single value from the subContBlk.
        
        deltaGLM = fit.GLMmulti(subTestGalvoBlk, 'simpLogSplitVSplitA');
        deltaGLM.prmInit = contParams{j}(i,:);
        deltaGLM.blockData.freeP = freeP;
        deltaGLM.fit;
        deltaParams{j}(i,:) = deltaGLM.prmFits;
        disp([i j]);
    end
end
savePath = which('spatialAnalysis');
savePath = savePath(1:strfind(savePath, '\@sp'));
savePath = [savePath '\data4Plots\ModelPerturbations\ModelPerturbationShuffles' datestr(now,'YYMMDD') '.mat'];

prmLabels = {'Bias'; 'visScaleIpsi'; 'visScaleConta'; 'N'; 'audScaleIpsi'; 'audScaleContra'};
save(savePath, 'normEstRepeats', 'contParams', 'deltaParams', 'gridXY', 'prmLabels', 'freeP');

inactResultsForModel.gridXY = gridXY;
inactResultsForModel.normEstRepeats = normEstRepeats;
inactResultsForModel.contParams = contParams;
inactResultsForModel.deltaParams = deltaParams;
inactResultsForModel.prmLabels = prmLabels;
inactResultsForModel.freeP = freeP;

%%
%Make the "scanPlot" structure for plotting. Note that "contData" is a constant here, so just subtract from the mean over hte inactivaiton
%grids for the different subsamples
%Set up plotting axes arrangement on figrure


axesOpt.btlrMargins =  [10 30 10 10];
axesOpt.gapBetweenAxes = [10 0];
axesOpt.axesSize = [200 200];

axesOpt.numOfRows = 2;
axesOpt.numOfCols = 3;
axesOpt.totalNumOfAxes = 6;

scanPlot.gridXY = gridXY;
for i = find(freeP)
    nShuffles = size(deltaParams{1},1) - normEstRepeats;
    contParams(cellfun(@isempty, contParams)) = deal({nan*ones(max(max(cellfun(@length, contParams))),6)});
    contData = cellfun(@(x) nanmean(x(1:normEstRepeats,i)),deltaParams);
    stdData = cellfun(@(x) nanstd(x(normEstRepeats+1:end,i)),deltaParams);
    sortedData = arrayfun(@(x,y) sort(abs([x{1}(normEstRepeats+1:end,i);y]),'descend'), deltaParams,contData, 'uni', 0);
    
    obj.hand.axes = plt.getAxes(axesOpt, i);
    scanPlot.title = prmLabels{i};
    scanPlot.data = contData./stdData;
    scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(contData), sortedData,'uni', 0));
    scanPlot.pVals
    scanPlot.colorBarLimits = [-30 30];
    sigLevels = (10.^(-2:-1:-10))';
    lastSigLevel = find(sigLevels>min(scanPlot.pVals(:)),1,'last');
    scanPlot.sigLevels = [0.01; 0.001; 0.0001];%sigLevels(max([1 lastSigLevel-2]):lastSigLevel);
    plt.scanningBrainEffects(scanPlot);
end
end
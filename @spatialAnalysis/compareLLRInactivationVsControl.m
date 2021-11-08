function inactLLRResults = compareLLRInactivationVsControl(obj, nShuffles, crossVal, contOnly)
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
if ~exist('nShuffles', 'var'); nShuffles = 0; end
if ~exist('crossVal', 'var'); crossVal = 2; end
if ~exist('contOnly', 'var'); contOnly = 0; end
freeP = [1 1 1 0 1 1]>0;

%Create "iBlk" (initialBlock) which removes some incorrect galvoPositions, repeated trials, and keeps only valid trials
iBlk = prc.filtBlock(obj.blks, obj.blks.tri.inactivation.galvoPosition(:,2)~=4.5);
iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));

%Conditional to optionally group inactivaiton sites together (e.g. if you want to combine all V1 sites). We "reflect" these groups, so only one
%hemisphere needs to be defined. We find all trials with galvoPositions belonging to those groups in iBlk and replace with the mean group position
posNames = {'MOs'; 'V1'; 'A1'; 'S1'};
galvoGrps = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4];[3,1; 3,0; 4.2,0]};
tstDat = [abs(iBlk.tri.inactivation.galvoPosition(:,1)) iBlk.tri.inactivation.galvoPosition(:,2)];
grpNum = length(galvoGrps);
grpIdx = cellfun(@(x,y) ismember(tstDat, x, 'rows').*y, galvoGrps, num2cell(1:grpNum)', 'uni', 0);
grpIdx = sum(cell2mat(grpIdx'),2);
meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
iBlk.tri.inactivation.galvoPosition(grpIdx>0,:) = cell2mat(meanPositions(grpIdx(grpIdx>0)));
iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | grpIdx>0);


contBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 0);
uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 1);

%If plotType contains 'dif' then we want to switch the responseCalc, vis and aud paramters such that all "left" trials are flipped (vis left trials in
%the case of conflict trials). Now, inactivations on the right hemisphere are contralateral, and left hemisphere is ipsilateral
idx2Flip = uniBlk.tri.inactivation.galvoPosition(:,1)<0;
uniBlk.tri.outcome.responseCalc(idx2Flip) = (uniBlk.tri.outcome.responseCalc(idx2Flip)*-1+3).*(uniBlk.tri.outcome.responseCalc(idx2Flip)>0);
uniBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*uniBlk.tri.inactivation.galvoPosition(idx2Flip,1);
uniBlk.tri.stim.audDiff(idx2Flip) = -1*uniBlk.tri.stim.audDiff(idx2Flip);
uniBlk.tri.stim.visDiff(idx2Flip) = -1*uniBlk.tri.stim.visDiff(idx2Flip);
uniBlk.tri.stim.conditionLabel(idx2Flip) = -1*uniBlk.tri.stim.conditionLabel(idx2Flip);


[~, gridXY] = prc.makeGrid(uniBlk, uniBlk.tri.outcome.responseCalc, [], 'galvouni',2);

%Define number of suffles, and the number of times to estimate the control (default is 10% or 500) since this will change with different
%subsamples of the subjetcs. Total loops is the sum of shuffle and control loops. We randomly assign galvoPositions to the contBlk from the
%galvoPositions in the testBlock (for shuffling purposes)
if ~nShuffles; nShuffles = 1500; end
nShuffles = round(nShuffles/10)*10;
if contOnly
    normEstRepeats = nShuffles;
else
    normEstRepeats = round(nShuffles/10);
end

%Use some matrix tricks to create "uniformLaserFilters" which is a filter for each shuffle that equalizes the frequency of subjects
%contributing to each point in the grid of galvo positions.
trialIdx = (1:uniBlk.tot.trials)';
nonUniformLaser = prc.makeGrid(uniBlk, [double(uniBlk.tri.subjectRef) trialIdx], [], 'galvouni',2);
laserShuffles = cellfun(@(x) double(prc.makeFreqUniform(x(:,1),normEstRepeats,x(:,2))), nonUniformLaser, 'uni', 0);
laserShuffles = num2cell(cell2mat(laserShuffles(:)),1)';
uniformLaserFilters = repmat(cellfun(@(x) ismember(trialIdx, x), laserShuffles, 'uni', 0),11,1);
uniformControlFilters = repmat(num2cell(prc.makeFreqUniform(contBlk.tri.subjectRef,normEstRepeats),1)',11,1);

%Removing these excess fields makes "filtBlock" run significantly faster in the subsequent loop
uniBlk.tri = rmfield(uniBlk.tri, {'timings'});
contBlk.tri = rmfield(contBlk.tri, {'timings'});
uniBlk.tri.inactivation = rmfield(uniBlk.tri.inactivation, {'laserPower', 'galvoType', 'laserOnsetDelay', 'laserDuration'});
contBlk.tri.inactivation = rmfield(contBlk.tri.inactivation, {'laserPower', 'galvoType', 'laserOnsetDelay', 'laserDuration'});

%This is the same loop as above, to generate "inactiveGrid" but with some extra steps to deal with the shuffles
[contParams, logLikDelta, logLikNoDelta, deltaParams] = deal(cell(length(posNames),1));

if contOnly; totalLoops = normEstRepeats; else, totalLoops = normEstRepeats+nShuffles; end
for i = 1:totalLoops
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
    
    for j = 4:length(posNames)
        contParams{j}(i,:) = contGLM.prmFits;
        galvoRef = meanPositions{j};
        includeIdx = ismember(subTestBlk.tri.inactivation.galvoPosition,galvoRef, 'rows');
        subTestGalvoBlk = prc.filtBlock(subTestBlk, includeIdx);
        
        %After shuffling, separate randomBlk based on laser activation. This will do nothing if unshuffled of course. Then estimate the result
        %(which is based on data2Use and op2Use) for each galvoPosition for subTestBlk and subtract the single value from the subContBlk.
        
        deltaGLM = fit.GLMmulti(subTestGalvoBlk, 'simpLogSplitVSplitA');
        deltaGLM.prmInit = contParams{j}(i,:);
        deltaGLM.blockData.freeP = freeP;
        deltaGLM.fitCV(crossVal);
        deltaParams{j}(i,:) = mean(deltaGLM.prmFits,1);
        logLikDelta{j}(i,:) = mean(deltaGLM.logLik);
        
        deltaGLM = fit.GLMmulti(subTestGalvoBlk, 'simpLogSplitVSplitA');
        deltaGLM.prmInit = contParams{j}(i,:);
        deltaGLM.blockData.freeP = (freeP*0)>1;
        deltaGLM.fitCV(crossVal);
        logLikNoDelta{j}(i,:) = mean(deltaGLM.logLik);

        disp([i j]);
    end
end
% savePath = 'D:\Dropbox (Neuropixels)\MouseData';
% savePath = [savePath '\LLRInativationSignificance' datestr(now,'YYMMDD') '.mat'];
% % 
prmLabels = {'Bias'; 'visScaleIpsi'; 'visScaleConta'; 'N'; 'audScaleIpsi'; 'audScaleContra'};
% if ~contOnly
%     save(savePath, 'normEstRepeats', 'contParams', 'deltaParams', 'gridXY', 'prmLabels', 'freeP','logLikDelta','logLikNoDelta');
% end

inactLLRResults.posNames = posNames;
inactLLRResults.meanPositions = meanPositions;
inactLLRResults.normEstRepeats = normEstRepeats;
inactLLRResults.contParams = contParams;
inactLLRResults.deltaParams = deltaParams;
inactLLRResults.logLikDelta = logLikDelta;
inactLLRResults.logLikNoDelta = logLikNoDelta;
inactLLRResults.prmLabels = prmLabels;
inactLLRResults.freeP = freeP;
end
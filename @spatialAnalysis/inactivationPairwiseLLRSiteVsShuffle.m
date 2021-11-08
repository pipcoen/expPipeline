function inactCompResults = inactivationPairwiseLLRSiteVsShuffle(obj, nShuffles, pInclude)
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
if ~exist('pInclude', 'var'); pInclude = 1:5; end
freeP = [1 1 1 1 1 1]>0;

%Create "iBlk" (initialBlock) which removes some incorrect galvoPositions, repeated trials, and keeps only valid trials
iBlk = prc.filtBlock(obj.blks, obj.blks.tri.inactivation.galvoPosition(:,2)~=4.5);
iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
iBlk.tri.inactivation.galvoPosition(iBlk.tri.inactivation.laserType==0,:) = repmat([0 0],sum(iBlk.tri.inactivation.laserType==0),1);

%Conditional to optionally group inactivaiton sites together (e.g. if you want to combine all V1 sites). We "reflect" these groups, so only one
%hemisphere needs to be defined. We find all trials with galvoPositions belonging to those groups in iBlk and replace with the mean group position
posNames = {'MOs'; 'V1'; 'A1'; 'S1';'None'};
prmLabels = {'Bias'; 'visScaleIpsi'; 'visScaleConta'; 'N'; 'audScaleIpsi'; 'audScaleContra'};
posNames = posNames(pInclude);
galvoGrps = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4];[3,1; 3,0; 4.2,0];[0 0 ; 0 0]};
galvoGrps = galvoGrps(pInclude);

tstDat = [abs(iBlk.tri.inactivation.galvoPosition(:,1)) iBlk.tri.inactivation.galvoPosition(:,2)];
grpNum = length(galvoGrps);
grpIdx = cellfun(@(x,y) ismember(tstDat, x, 'rows').*y, galvoGrps, num2cell(1:grpNum)', 'uni', 0);
grpIdx = sum(cell2mat(grpIdx'),2);
meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
iBlk.tri.inactivation.grpIdx = grpIdx;
iBlk.tri.inactivation.galvoPosition(grpIdx>0,:) = cell2mat(meanPositions(grpIdx(grpIdx>0)));
iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | grpIdx>0);

nBlk = iBlk;

%If plotType contains 'dif' then we want to switch the responseCalc, vis and aud paramters such that all "left" trials are flipped (vis left trials in
%the case of conflict trials). Now, inactivations on the right hemisphere are contralateral, and left hemisphere is ipsilateral
idx2Flip = nBlk.tri.inactivation.galvoPosition(:,1)<0;
nBlk.tri.outcome.responseCalc(idx2Flip) = (nBlk.tri.outcome.responseCalc(idx2Flip)*-1+3).*(nBlk.tri.outcome.responseCalc(idx2Flip)>0);
nBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*nBlk.tri.inactivation.galvoPosition(idx2Flip,1);
nBlk.tri.stim.audDiff(idx2Flip) = -1*nBlk.tri.stim.audDiff(idx2Flip);
nBlk.tri.stim.visDiff(idx2Flip) = -1*nBlk.tri.stim.visDiff(idx2Flip);
nBlk.tri.stim.conditionLabel(idx2Flip) = -1*nBlk.tri.stim.conditionLabel(idx2Flip);

%Define number of suffles, and the number of times to estimate the control (default is 10% or 500) since this will change with different
%subsamples of the subjetcs. Total loops is the sum of shuffle and control loops. We randomly assign galvoPositions to the contBlk from the
%galvoPositions in the testBlock (for shuffling purposes)
if ~nShuffles; nShuffles = 1500; end
nShuffles = round(nShuffles/10)*10;
normEstRepeats = round(nShuffles/10);

%Use some matrix tricks to create "uniformLaserFilters" which is a filter for each shuffle that equalizes the frequency of subjects
%contributing to each point in the grid of galvo positions.
trialIdx = (1:nBlk.tot.trials)';
nonUniformLaser = prc.makeGrid(nBlk, [double(nBlk.tri.subjectRef) trialIdx], [], 'galvouni',2);
laserShuffles = cellfun(@(x) double(prc.makeFreqUniform(x(:,1),normEstRepeats,x(:,2))), nonUniformLaser, 'uni', 0);
laserShuffles = num2cell(cell2mat(laserShuffles(:)),1)';
uniformLaserFilters = repmat(cellfun(@(x) ismember(trialIdx, x), laserShuffles, 'uni', 0),11,1);

%Removing these excess fields makes "filtBlock" run significantly faster in the subsequent loop
nBlk.tri = rmfield(nBlk.tri, {'timings'});
nBlk.tri.inactivation = rmfield(nBlk.tri.inactivation, {'laserPower', 'galvoType', 'laserOnsetDelay', 'laserDuration'});

%This is the same loop as above, to generate "inactiveGrid" but with some extra steps to deal with the shuffles

totalLoops = normEstRepeats+nShuffles;
[grpIdx1, grpIdx2] = meshgrid(1:5, 1:5);
grpOrd = [grpIdx1(:) grpIdx2(:)];
grpOrd(diff(grpOrd,[],2)==0,:) = [];

tDat = nBlk.tri.inactivation.galvoPosition;
mPos = cell2mat(meanPositions);
grpFilter = cellfun(@(x) ismember(tDat, mPos([x(1) x(2)],:), 'rows'), num2cell(grpOrd,2), 'uni', 0);

[trainParams, testParams, logLikTest, logLikTrain, selfTestParams, logLikSelfTest] = deal(cell(size(grpOrd,1),1));
for j = 1:size(grpOrd,1)
    for i = 1:totalLoops
        %Filter both blks to make contributions from each subject equal at each galvoPosition and across control trials.
        uniformBlk = prc.filtBlock(nBlk, uniformLaserFilters{i} & grpFilter{j});
        
        %Generate "randomBlk" by concatenating control and test blocks. If the repeat number is outside the "normEstRepeats" range then we shuffle
        %the laser activation and galvoPositions
        if i > normEstRepeats
            uniformBlk.tri.inactivation.grpIdx = uniformBlk.tri.inactivation.grpIdx(randperm(uniformBlk.tot.trials),:);
        end
        trainBlk = prc.filtBlock(uniformBlk, uniformBlk.tri.inactivation.grpIdx == grpOrd(j,1));
        testBlk = prc.filtBlock(uniformBlk, uniformBlk.tri.inactivation.grpIdx == grpOrd(j,2));
        
        trainGLM = fit.GLMmulti(trainBlk, 'simpLogSplitVSplitA');
        testGLM.blockData.freeP = freeP;
        trainGLM.fit
        trainParams{j}(i,:) = trainGLM.prmFits;
        logLikTrain{j}(i,:) = mean(trainGLM.logLik);
        
        %After shuffling, separate randomBlk based on laser activation. This will do nothing if unshuffled of course. Then estimate the result
        %(which is based on data2Use and op2Use) for each galvoPosition for testBlk and subtract the single value from the trainBlk.        
        testGLM = fit.GLMmulti(testBlk, 'simpLogSplitVSplitA');
        testGLM.prmInit = trainParams{j}(i,:);
        testGLM.blockData.freeP = (freeP*0)>0;
        testGLM.fit;
        testParams{j}(i,:) = mean(testGLM.prmFits,1);
        logLikTest{j}(i,:) = mean(testGLM.logLik);
        
        selfTestGLM = fit.GLMmulti(testBlk, 'simpLogSplitVSplitA');
        selfTestGLM.prmInit = trainParams{j}(i,:);
        selfTestGLM.blockData.freeP = freeP;
        selfTestGLM.fit;
        selfTestParams{j}(i,:) = mean(selfTestGLM.prmFits,1);
        logLikSelfTest{j}(i,:) = mean(selfTestGLM.logLik);
                
        disp([i j]);
    end
end
%%
inactCompResults.posNames = posNames;
inactCompResults.grpOrd = grpOrd;
inactCompResults.prmLabels = prmLabels;
inactCompResults.meanPositions = meanPositions;
inactCompResults.normEstRepeats = normEstRepeats;
inactCompResults.trainParams = trainParams;
inactCompResults.testParams = testParams;
inactCompResults.selfTestParams = selfTestParams;
inactCompResults.logLikTest = logLikTest;
inactCompResults.logLikTrain = logLikTrain;
inactCompResults.logLikSelfTest = logLikSelfTest;

%%
savePath = 'D:\Dropbox (Neuropixels)\MouseData';
savePath = [savePath '\inactCompResults' datestr(now,'YYMMDD') '.mat'];
save(savePath, '-struct', 'inactCompResults')
end
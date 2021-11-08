function [logLik, fitRegion, testRegion] = compareLLRInactivationAcrossRegions(obj)
%% Method for "spatialAnalysis" class. Plots effects of inactivation on behavior. Plots are shown as grids on an outline of cortex.

%INPUTS(default values)
%plotType(res)---------A string that contains three letter tags indicates the type of plot to generate. Can be combination of below options
%   'res'-------------------------quantify changes in the mouse response (i.e. fraciton of rightward choices)
%   'dif'-------------------------rather than separately analysing left and right trials, combine trials and use ipsilateral and contralateral
%   'grp'-------------------------combine inactivation sites into M2 and Vis
%   'sig'-------------------------test the significance of inactivation by shuffling inactivation sites and laser on/off
%nRepeats-------------The number of times to shuffle te data
%subsets---------------The data subsets to make plots for (see function "prc.getDefinedSubset")

%Set up defaults for the input values. "op2use" is "mean" as default for responses, but changes below depending on the type data being used
numOfMice = length(obj.blks);
if numOfMice > 1; error('Only coded to handle one mouse atm'); end
freeP = [1 1 1 1 1 1]>0;

%Create "iBlk" (initialBlock) which removes some incorrect galvoPositions, repeated trials, and keeps only valid trials
iBlk = prc.filtBlock(obj.blks, obj.blks.tri.inactivation.galvoPosition(:,2)~=4.5);
iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));

%If plotType contains 'dif' then we want to switch the responseCalc, vis and aud paramters such that all "left" trials are flipped (vis left trials in
%the case of conflict trials). Now, inactivations on the right hemisphere are contralateral, and left hemisphere is ipsilateral
idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType == 1;
iBlk.tri.outcome.responseCalc(idx2Flip) = (iBlk.tri.outcome.responseCalc(idx2Flip)*-1+3).*(iBlk.tri.outcome.responseCalc(idx2Flip)>0);
iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);

%Conditional to optionally group inactivaiton sites together (e.g. if you want to combine all V1 sites). We "reflect" these groups, so only one
%hemisphere needs to be defined. We find all trials with galvoPositions belonging to those groups in iBlk and replace with the mean group position
posNames = {'MOs'; 'V1'; 'A1'; 'S1'; 'None'};
galvoGrps = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4];[3,1; 3,0; 4.2,0]};
tstDat = [abs(iBlk.tri.inactivation.galvoPosition(:,1)) iBlk.tri.inactivation.galvoPosition(:,2)];
grpIdx = cellfun(@(x,y) ismember(tstDat, x, 'rows').*y, galvoGrps, num2cell(1:length(galvoGrps))', 'uni', 0);
grpIdx = sum(cell2mat(grpIdx'),2);
meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
iBlk.tri.inactivation.galvoPosition(grpIdx>0,:) = cell2mat(meanPositions(grpIdx(grpIdx>0)));
iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | grpIdx>0);

contBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 0);
uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 1);

%%
[fitRegion, testRegion] = deal(cell(5,5)); 
logLik = nan*ones(5,5);
for i = 1:5
    for j = 1:5
        if i == 5
            trainBlk = contBlk;
        else
            trainBlk = prc.filtBlock(uniBlk, ismember(abs(uniBlk.tri.inactivation.galvoPosition), abs(meanPositions{i}), 'rows'));
        end
        fitRegion(i,j) = posNames(i);
        if j == 5
            testBlk = contBlk;
        else
            testBlk = prc.filtBlock(uniBlk, ismember(abs(uniBlk.tri.inactivation.galvoPosition), abs(meanPositions{j}), 'rows'));
        end
        testRegion(i,j) = posNames(j);
        if trainBlk.tot.trials > floor(0.5*testBlk.tot.trials)
            tDat = zeros(trainBlk.tot.trials,1) > 0;
            tDat(randperm(trainBlk.tot.trials, floor(testBlk.tot.trials/2))) = 1;
            trainBlk = prc.filtBlock(trainBlk, tDat);
        else
            tDat = zeros(testBlk.tot.trials,1) > 0;
            tDat(randperm(testBlk.tot.trials, floor(trainBlk.tot.trials*2))) = 1;
            testBlk = prc.filtBlock(testBlk, tDat);
        end
        
        
        trainGLM = fit.GLMmulti(trainBlk, 'simpLogSplitVSplitA');
        trainGLM.blockData.freeP = freeP;
        trainGLM.fit;
        
        testGLM = fit.GLMmulti(testBlk, 'simpLogSplitVSplitA');
        testGLM.prmInit = mean(trainGLM.prmFits,1);
        testGLM.blockData.freeP = (freeP*0)>0;
        testGLM.fitCV(2);
        logLik(i,j) = mean(testGLM.logLik);

        disp([i j]);
    end
end
inactLLRResults.logLik = logLik;
inactLLRResults.fitRegion = fitRegion;
inactLLRResults.testRegion = testRegion;
end
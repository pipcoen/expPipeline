function visAndM2Inactivation
%% This function plots the data panels for figure one of the ms
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

uniData = spatialAnalysis('all', 'uniscan', 1, 1);
iBlk = prc.filtBlock(uniData.blks, uniData.blks.tri.inactivation.galvoPosition(:,2)~=4.5);
iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, iBlk.tri.outcome.responseMade~=0);

galvoGrps = {[1.8 -4; 3,-4; 3,-3];[0.6 2; 1.8, 2; 0.6, 3]};
galvoGrps = [galvoGrps; cellfun(@(x) [x(:,1)*-1, x(:,2)], galvoGrps, 'uni', 0)];
grpNum = length(galvoGrps);
grpIdx = cellfun(@(x,y) ismember(iBlk.tri.inactivation.galvoPosition, x, 'rows').*y, galvoGrps, num2cell(1:grpNum)', 'uni', 0);
grpIdx = sum(cell2mat(grpIdx'),2);
meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
iBlk.tri.inactivation.galvoPosition(grpIdx>0,:) = cell2mat(meanPositions(grpIdx(grpIdx>0)));
iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | grpIdx>0);

idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0;
iBlk.tri.outcome.responseMade(idx2Flip) = (iBlk.tri.outcome.responseMade(idx2Flip)*-1+3).*(iBlk.tri.outcome.responseMade(idx2Flip)>0);
iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
prmLabels = {'Bias'; 'visScaleIpsi'; 'visScaleConta'; 'N'; 'audScaleIpsi'; 'audScaleContra'};
%%
contBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 0);
nContShuffles = 1;
contParams = zeros(nContShuffles,6);
for i = 1:nContShuffles
    tempBlk = prc.filtBlock(contBlk, prc.makeFreqUniform(contBlk.tri.subjectRef));
    contGLM = fit.GLMmulti(tempBlk, 'simpLogSplitVSplitA');
    contGLM.fit;
    contParams(i,:) = contGLM.prmFits;
end
contGLM.prmFits = contParams(i,:);

contObj = copy(uniData);
contObj.blks = contBlk;
contObj.glmFit{1} = contGLM;


laserBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType == 1);
galvoRef = unique(laserBlk.tri.inactivation.galvoPosition, 'rows');
%%
figure;
freeP = [1 0 0 0 0 0; 0 0 1 0 0 0];
for i = 1:size(galvoRef,1)
    axesHandle = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
    tempBlk = prc.filtBlock(laserBlk, ismember(laserBlk.tri.inactivation.galvoPosition,galvoRef(i,:), 'rows'));
    lasObj = copy(uniData);
    lasObj.blks = tempBlk;
    lasObj.glmFit{1} = contGLM;
    lasObj.viewGLMFits('simpLogSplitVSplitA', [], 'only', 1)
    axis square;
    
    lasObj.blks.freeP = freeP(i,:)>0;
    lasObj.glmFit{1} = fit.GLMmulti(lasObj.blks, 'simpLogSplitVSplitA');
    lasObj.glmFit{1}.prmInit = contGLM.prmFits;
    lasObj.glmFit{1}.fit;
%     lasObj.glmFit{1}.prmFits = lasObj.glmFit{1}.prmFits.*(freeP(i,:)>0)+lasObj.glmFit{1}.prmInit;
    axesHandle = plt.tightSubplot(nRows,nCols,i+size(galvoRef,1),axesGap,botTopMarg,lftRgtMarg);
    lasObj.viewGLMFits('simpLogSplitVSplitA', [], 'only', 1)
    axis square;
end

%%
export_fig('D:\Dropbox (Personal)\TalksAndApps\Papers\Coen_2020\3_visAndM2Inactivation', '-pdf', '-painters');
end
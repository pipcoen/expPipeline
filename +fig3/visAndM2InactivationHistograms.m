function visAndM2InactivationHistograms(groups)
%% This function plots the data panels for figure one of the ms
s = spatialAnalysis('all', 'uniscan', 1, 1);
if ~exist('groups', 'var'); groups = 'mosa1v1'; end
if ~exist('subsets', 'var'); subsets = {'VL', 'AL', 'CohL', 'ConL'}; end

figure;
axHeight = 250;
axWidth = 250;
nCols = 4;
nRows = 2;
figHeight = nRows*axHeight;
figWidth = nCols*axWidth;

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 100 figWidth, figHeight]);

%%
iBlk = prc.filtBlock(s.blks, s.blks.tri.inactivation.galvoPosition(:,2)~=4.5 | s.blks.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial);


groups = lower(groups);
[galvoGrps, grpName, grpColor] = deal({});
if contains(groups, 'v1')
    galvoGrps = [galvoGrps;  [1.8 -4; 3,-4; 3,-3]]; 
    grpName = [grpName; 'V1']; 
    grpColor = [grpColor; [0.9290, 0.6940, 0.1250]	];
end
if contains(groups, 'a1')
    galvoGrps = [galvoGrps; [4.2,-2; 4.2,-3; 4.2,-4]];
    grpName = [grpName; 'A1']; 
    grpColor = [grpColor; 'm'];
end
if contains(groups, 'mos')
    galvoGrps = [galvoGrps; [0.6 2; 1.8, 2; 0.6, 3]]; 
    grpName = [grpName; 'MOs']; 
    grpColor = [grpColor; [0.9100    0.4100    0.1700]];
end
if contains(groups, 'av')
    galvoGrps = [galvoGrps; [galvoGrps; [1.8 -4; 3,-4; 4.2,-4; 1.8,-3; 3,-3; 4.2,-3; 1.8,-2; 3,-2; 4.2,-2]]]; 
    grpName = [grpName; 'AV']; 
    grpColor = [grpColor; 'm'];
end

galvoGrps = [galvoGrps; cellfun(@(x) [x(:,1)*-1, x(:,2)], galvoGrps, 'uni', 0)];
grpNum = length(galvoGrps);
grpIdx = cellfun(@(x,y) ismember(iBlk.tri.inactivation.galvoPosition, x, 'rows').*y, galvoGrps, num2cell(1:grpNum)', 'uni', 0);
grpIdx = sum(cell2mat(grpIdx'),2);
meanPositions = cellfun(@mean, galvoGrps, 'uni', 0);
iBlk.tri.inactivation.galvoPosition(grpIdx>0,:) = cell2mat(meanPositions(grpIdx(grpIdx>0)));
iBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0 | grpIdx>0);

rIdx = iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0);
iBlk.tri.outcome.responseMade(rIdx) = (iBlk.tri.outcome.responseMade(rIdx)*-1+3).*(iBlk.tri.outcome.responseMade(rIdx)>0);
iBlk.tri.stim.audInitialAzimuth(rIdx) = iBlk.tri.stim.audInitialAzimuth(rIdx)*-1;
iBlk.tri.stim.visInitialAzimuth(rIdx) = iBlk.tri.stim.visInitialAzimuth(rIdx)*-1;
iBlk.tri.stim.visInitialAzimuth(isinf(iBlk.tri.stim.visInitialAzimuth)) = inf;
iBlk.tri.inactivation.galvoPosition(rIdx,1) = -1*iBlk.tri.inactivation.galvoPosition(rIdx,1);
iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.threshMoveTime));

normBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0);
uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==1);
%%
titles = {'Visual Contra'; 'Auditory  Contra'; 'Coherent Contra'; 'Conflict Contra'; ...
    'Visual Ipsi'; 'Auditory  Ipsi'; 'Coherent Ipsi'; 'Conflict Ipsi'};
nSubs = length(subsets);
for i  = 1:nSubs*2
    %Get subset of trials based on tag. e.g. 'vl' for visual left for both normBlk (contBlk) and uniBlk (testBlk). We need these irrespective of 'sig'
    pIdx = mod(i,4) + (mod(i,4)==0)*4;
    contBlk = prc.getDefinedSubset(normBlk, subsets{pIdx});
    contBlk = prc.filtBlock(contBlk, prc.makeFreqUniform(contBlk.tri.subjectRef));
    testBlk = prc.getDefinedSubset(uniBlk, subsets{pIdx});
    testBlk = prc.filtBlock(testBlk, prc.makeFreqUniform(testBlk.tri.subjectRef));
    plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
    xLimits = [100 400];
    binEdges = xLimits(1):10:xLimits(2);
    xDat = xLimits(1)+5:10:xLimits(2);
    
    contFrq = histcounts(contBlk.tri.outcome.threshMoveTime*1000, binEdges, 'normalization', 'cdf');
    plot(xDat, contFrq, 'k', 'linewidth', 2);
    for j = 1:length(grpName)
        if i > nSubs; continue; end
        hold on;
        subTestBlk = prc.filtBlock(testBlk, ismember(testBlk.tri.inactivation.galvoPosition, meanPositions{j}, 'rows'));
        tstFreq = histcounts(subTestBlk.tri.outcome.threshMoveTime*1000, binEdges, 'normalization', 'cdf');
        plot(xDat, tstFreq, 'color', grpColor{j}, 'linewidth', 2);
    end
   
    for j = length(grpName)+1:length(grpName)*2
        hold on;
        if i < nSubs+1; continue; end
        subTestBlk = prc.filtBlock(testBlk, ismember(testBlk.tri.inactivation.galvoPosition, meanPositions{j}, 'rows'));
        tstFreq = histcounts(subTestBlk.tri.outcome.threshMoveTime*1000, binEdges, 'normalization', 'cdf');
        plot(xDat, tstFreq, 'color', grpColor{j-3}, 'linewidth', 2);
    end
    plot(xlim, [0.5 0.5], '--k')
    box off;
    set(gca, 'xTick', xLimits)
    set(gca, 'yTick', [0 1])
    title(titles{i});
%     legend({'Cont', 'Vis', 'Aud', 'MOs'})
end
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_modelInactivationReactionTimeCDFS', '-pdf', '-painters');

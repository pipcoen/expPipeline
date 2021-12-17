function viewInactivationTimedEffects(obj, inactvationSite, trialType, plotType, contra, yRng)
%% Method for "spatialAnalysis" class. Plots effects of inactivation on behavior. Plots are shown as grids on an outline of cortex.

%INPUTS(default values)
%plotType(res)---------A string that contains three letter tags indicates the type of plot to generate. Can be combination of below options
%   'res'-------------------------quantify changes in the mouse response (i.e. fraciton of rightward choices)
%   'dif'-------------------------rather than separately analysing left and right trials, combine trials and use ipsilateral and contralateral
%   'grp'-------------------------combine inactivation inactvationSite into M2 and Vis
%   'sig'-------------------------test the significance of inactivation by shuffling inactivation inactvationSite and laser on/off
%nShuffles-------------The number of times to shuffle te data
%subsets---------------The data subsets to make plots for (see function "prc.getDefinedSubset")

%Set up defaults for the input values. "op2use" is "mean" as default for responses, but changes below depending on the type data being used
if exist('inactvationSite', 'var') && isempty(inactvationSite); clear inactvationSite; end
if ~exist('inactvationSite', 'var'); inactvationSite = 'mos'; end
if ~exist('contra', 'var'); contra = 1; end
if ~exist('trialType', 'var'); trialType = 'vis'; end
if ~exist('plotType', 'var'); plotType = 'res'; end
if ~exist('yRng', 'var'); yRng = [-0.1 0.4]; end


%%
binSize = 70;
plotRange = [-100-binSize/2 150+binSize/2];

overlap = binSize-10;
bins = (plotRange(1):(binSize-overlap):plotRange(2)-binSize)';
bins = [bins bins+binSize];
pLim = 10^-3;
xDat = bins(:,1)+binSize/2;

regRef = {'mos', [0.5, 2.0]; 'v1', [2.0, -4.0]; 's1' , [3.5 -0.5]; 'out' , [3 5.5]};

%Set up plotting axes arrangement on figrure
inactvationSite2keep = cell2mat(regRef(contains(regRef(:,1), inactvationSite),2));

iBlk = prc.filtBlock(obj.blks, obj.blks.tri.trialType.repeatNum==1 & obj.blks.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, ismember(abs(iBlk.tri.inactivation.galvoPosition),abs(inactvationSite2keep), 'rows') | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.stim.visContrast > 0 | iBlk.tri.trialType.auditory);

fIdx = iBlk.tri.trialType.visual*0;
fIdx = fIdx + contains(trialType, 'vis')*iBlk.tri.trialType.visual ...
    + contains(trialType, 'aud')*iBlk.tri.trialType.auditory ...
    + contains(trialType, 'coh')*iBlk.tri.trialType.coherent ...
    + contains(trialType, 'con')*iBlk.tri.trialType.conflict;
iBlk = prc.filtBlock(iBlk, fIdx);
if strcmpi(trialType, 'vis');  pltCol = [0.9290, 0.6940, 0.1250]; end
if strcmpi(trialType, 'aud');  pltCol = [189, 48, 255]/255; end

idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType==1;
iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
iBlk.tri.outcome.responseCalc(idx2Flip) = (iBlk.tri.outcome.responseCalc(idx2Flip)*-1+3).*(iBlk.tri.outcome.responseCalc(idx2Flip)>0);

if contains(plotType, 'res')
    iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
    iBlk.tri.inactivation.data2Use = iBlk.tri.outcome.responseCalc==2;
    vNames = {'right','left'};
elseif contains(plotType, 'tim')
    iBlk.tri.inactivation.data2Use = iBlk.tri.outcome.responseRecorded==0;
    vNames = {'timeout','response'};
elseif contains(plotType, 'rea')
    iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
    iBlk.tri.inactivation.data2Use = iBlk.tri.outcome.reactionTime;    
end
iBlk = prc.filtBlock(iBlk, iBlk.exp.numOfTrials>75);

if contra
    rIdx = iBlk.tri.stim.visDiff<0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff<0);
else
    rIdx = iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0);
end
iBlk = prc.filtBlock(iBlk, rIdx);

%Create normBlk and uniBlk which are filtered versions of iBlk with only control or inactivation trials respectively
normBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0);
uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==1);
laserOffsets = uniBlk.tri.inactivation.laserOnsetDelay*1000;
sampIdx = arrayfun(@(x,y) laserOffsets>x & laserOffsets<y, bins(:,1), bins(:,2), 'uni', 0);

if ~contains(plotType, 'rea')
    nrmDat.fracT = mean(normBlk.tri.inactivation.data2Use);
    nrmDat.nTrue = sum(normBlk.tri.inactivation.data2Use);
    nrmDat.nFalse = sum(~normBlk.tri.inactivation.data2Use);
    
    lasDat.fracT = cellfun(@(x) mean(uniBlk.tri.inactivation.data2Use(x)), sampIdx);
    lasDat.nTrue = cellfun(@(x) sum(uniBlk.tri.inactivation.data2Use(x)), sampIdx);
    lasDat.nFalse = cellfun(@(x) sum(~uniBlk.tri.inactivation.data2Use(x)), sampIdx);
    
    CI = 1.96*sqrt((lasDat.fracT.*(1-lasDat.fracT))./(lasDat.nTrue+lasDat.nFalse));
    pltM = lasDat.fracT-nrmDat.fracT;
else
    nrmDat.medR = median(normBlk.tri.inactivation.data2Use);
    lasDat.medR = cellfun(@(x) median(uniBlk.tri.inactivation.data2Use(x)), sampIdx);
    
    CI = cellfun(@(x) mad(uniBlk.tri.inactivation.data2Use(x)), sampIdx);
    pltM = lasDat.medR-nrmDat.medR;
end

plotData = permute(cat(3, pltM, pltM-CI, pltM+CI), [2 1 3]);
opt.Marker = 'none';
plt.rowsOfGrid(xDat', plotData, pltCol, opt);

if ~contains(plotType, 'rea')
    for i = 1:length(xDat)
        tbl = table([nrmDat.nTrue;lasDat.nTrue(i)],[nrmDat.nFalse;lasDat.nFalse(i)], ...
            'VariableNames',vNames,'RowNames',{'Cnt','Las'});
        [~,pVal(i)] = fishertest(tbl);
    end
else
    pVal = arrayfun(@(x) ranksum(x,normBlk.tri.inactivation.data2Use), lasDat.medR);
end
        
ylim(yRng)
plot(xlim, [0,0], '--k', 'linewidth', 1.5)
plot([0,0], ylim, '--k', 'linewidth', 1.5)
if any(pVal<pLim); plot(xDat(pVal<pLim), max(pltM+CI)+0.05, '.', 'color',  pltCol); end
end
function viewInactivationTimedEffectsRT(obj, sites, align, plotType)
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
if exist('sites', 'var') && isempty(sites); clear sites; end
if ~exist('sites', 'var'); sites = 'mos'; end
if ~exist('align', 'var'); align = 'move'; end
if ~exist('plotType', 'var'); plotType = 1; end
if strcmpi(align, 'stim'); stim = 1; else, stim = 0; end


%%
if stim
    plotRange = [-125 200];
    yRng = [-0.1 0.2];
else
    plotRange = [-300 100];
    yRng = [-0.1 0.15];
    if strcmpi(sites, 'v1'); yRng = [-0.3 0.45]; end
end

binSize = 50;
overlap = 48;
bins = (plotRange(1):(binSize-overlap):plotRange(2)-binSize)';
bins = [bins bins+binSize];
pLim = 10^-5;
xDat = bins(:,1);

regRef = {'mos', [0.5, 2.0]; 'v1', [2.0, -4.0]};

%Set up plotting axes arrangement on figrure
sites2keep = cell2mat(regRef(contains(regRef(:,1), sites),2));

iBlk = prc.filtBlock(obj.blks, obj.blks.tri.trialType.repeatNum==1 & obj.blks.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, ismember(abs(iBlk.tri.inactivation.galvoPosition),abs(sites2keep), 'rows') | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.stim.visContrast > 0 | iBlk.tri.trialType.auditory);
if strcmpi(sites, 'v1') 
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.visual); 
    pltCol = [0.9290, 0.6940, 0.1250];
end
if strcmpi(sites, 'mos')
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.visual);
    pltCol = [0.8500, 0.3250, 0.0980];
end

%If plotType contains 'dif' then we want to switch the responseCalc, vis and aud paramters such that all "right" trials are flipped (vis right trials in
%the case of conflict trials). Now, inactivations on the right hemisphere are contralateral, and left hemisphere is ipsilateral
rIdx = iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0);
iBlk.tri.outcome.responseCalc(rIdx) = (iBlk.tri.outcome.responseCalc(rIdx)*-1+3).*(iBlk.tri.outcome.responseCalc(rIdx)>0);
iBlk.tri.stim.audInitialAzimuth(rIdx) = iBlk.tri.stim.audInitialAzimuth(rIdx)*-1;
iBlk.tri.stim.visInitialAzimuth(rIdx) = iBlk.tri.stim.visInitialAzimuth(rIdx)*-1;
iBlk.tri.stim.visInitialAzimuth(isinf(iBlk.tri.stim.visInitialAzimuth)) = inf;
iBlk.tri.inactivation.galvoPosition(rIdx,1) = -1*iBlk.tri.inactivation.galvoPosition(rIdx,1);

iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime));
iBlk.tri.inactivation.data2Use = iBlk.tri.outcome.reactionTime*1000;

%Create normBlk and uniBlk which are filtered versions of iBlk with only control or inactivation trials respectively
normBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0);
uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==1);

nrmDat.fracR = mean(normBlk.tri.inactivation.data2Use);
nrmDat.nTri = length(normBlk.tri.inactivation.data2Use);
nrmDat.nRT = sum(normBlk.tri.inactivation.data2Use);

laserOffsets = uniBlk.tri.inactivation.laserOnsetDelay*1000;
if ~stim; laserOffsets = laserOffsets - uniBlk.tri.outcome.reactionTime*1000; end
sampIdx = arrayfun(@(x,y) laserOffsets>x & laserOffsets<y, bins(:,1), bins(:,2), 'uni', 0);
contra = uniBlk.tri.inactivation.galvoPosition(:,1)>0;

xlim([xDat(1) xDat(end)])
% ylim(yRng)
for i = 1:2
    hold on;
    if i == 2; contra = ~contra; end
    tDat.fracR(i,:) = cellfun(@(x) mean(uniBlk.tri.inactivation.data2Use(x & contra)), sampIdx);
    tDat.nTri(i,:) = cellfun(@(x) length(uniBlk.tri.inactivation.data2Use(x & contra)), sampIdx);
    tDat.nRT(i,:) = cellfun(@(x) sum(uniBlk.tri.inactivation.data2Use(x & contra)), sampIdx);
end

if plotType==1
    for i = 1:2
        pltM = tDat.fracR(i,:);
        plot(xDat', pltM, 'color', pltCol/i);
    end
end
end
function viewInactivationTimedEffects(obj, sites, align, triT, plotType, conT)
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
if ~exist('conT', 'var'); conT = 1; end
if ~exist('plotType', 'var'); plotType = 1; end
if ~exist('triT', 'var'); triT = 'vis'; end
if strcmpi(align, 'stim'); stim = 1; else, stim = 0; end


%%
binSize = 70;
if stim
    plotRange = [-100-binSize/2 150+binSize/2];
    if strcmpi(triT, 'vis'); yRng = [-0.1 0.4]; end
    if strcmpi(triT, 'aud'); yRng = [-0.1 0.4]; end
else
    plotRange = [-300 100];
    yRng = [-0.1 0.15];
    if strcmpi(sites, 'v1'); yRng = [-0.3 0.45]; end
end

overlap = binSize-10;
bins = (plotRange(1):(binSize-overlap):plotRange(2)-binSize)';
bins = [bins bins+binSize];
pLim = 10^-3;
xDat = bins(:,1)+binSize/2;

regRef = {'mos', [0.5, 2.0]; 'v1', [2.0, -4.0]; 's1' , [3.5 -0.5]; 'out' , [3 5.5]};

%Set up plotting axes arrangement on figrure
sites2keep = cell2mat(regRef(contains(regRef(:,1), sites),2));

iBlk = prc.filtBlock(obj.blks, obj.blks.tri.trialType.repeatNum==1 & obj.blks.tri.trialType.validTrial);
iBlk = prc.filtBlock(iBlk, ismember(abs(iBlk.tri.inactivation.galvoPosition),abs(sites2keep), 'rows') | iBlk.tri.inactivation.laserType==0);
iBlk = prc.filtBlock(iBlk, iBlk.tri.stim.visContrast > 0 | iBlk.tri.trialType.auditory);

fIdx = iBlk.tri.trialType.visual*0;
fIdx = fIdx + contains(triT, 'vis')*iBlk.tri.trialType.visual ...
    + contains(triT, 'aud')*iBlk.tri.trialType.auditory ...
    + contains(triT, 'coh')*iBlk.tri.trialType.coherent ...
    + contains(triT, 'con')*iBlk.tri.trialType.conflict;
iBlk = prc.filtBlock(iBlk, fIdx);
pltCol = [0 0 0];
if strcmpi(triT, 'vis') 
    pltCol = [0.9290, 0.6940, 0.1250];
end
if strcmpi(triT, 'aud')
    pltCol = [189, 48, 255]/255;
end

%If plotType contains 'dif' then we want to switch the responseCalc, vis and aud paramters such that all "right" trials are flipped (vis right trials in
%the case of conflict trials). Now, inactivations on the right hemisphere are contralateral, and left hemisphere is ipsilateral
% rIdx = iBlk.tri.stim.visDiff>0 | (iBlk.tri.stim.visDiff==0 & iBlk.tri.stim.audDiff>0);
% iBlk.tri.outcome.responseCalc(rIdx) = (iBlk.tri.outcome.responseCalc(rIdx)*-1+3).*(iBlk.tri.outcome.responseCalc(rIdx)>0);
% iBlk.tri.inactivation.galvoPosition(rIdx,1) = -1*iBlk.tri.inactivation.galvoPosition(rIdx,1);

idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType==1;
iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
iBlk.tri.outcome.responseCalc(idx2Flip) = (iBlk.tri.outcome.responseCalc(idx2Flip)*-1+3).*(iBlk.tri.outcome.responseCalc(idx2Flip)>0);

iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.responseCalc));
iBlk.tri.inactivation.data2Use = iBlk.tri.outcome.responseCalc==2;
iBlk = prc.filtBlock(iBlk, iBlk.exp.numOfTrials>75);


%Create normBlk and uniBlk which are filtered versions of iBlk with only control or inactivation trials respectively
normBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==0);
if conT
    rIdx = normBlk.tri.stim.visDiff<0 | (normBlk.tri.stim.visDiff==0 & normBlk.tri.stim.audDiff<0);
else
    rIdx = normBlk.tri.stim.visDiff>0 | (normBlk.tri.stim.visDiff==0 & normBlk.tri.stim.audDiff>0);
end

nrmDat.fracR = mean(normBlk.tri.inactivation.data2Use(rIdx));
nrmDat.nTri = length(normBlk.tri.inactivation.data2Use(rIdx));
nrmDat.nRT = sum(normBlk.tri.inactivation.data2Use(rIdx));

uniBlk = prc.filtBlock(iBlk, iBlk.tri.inactivation.laserType==1);
laserOffsets = uniBlk.tri.inactivation.laserOnsetDelay*1000;
if ~stim; laserOffsets = laserOffsets - uniBlk.tri.outcome.reactionTime*1000; end
sampIdx = arrayfun(@(x,y) laserOffsets>x & laserOffsets<y, bins(:,1), bins(:,2), 'uni', 0);

if conT
    contra = uniBlk.tri.stim.visDiff<0 | (uniBlk.tri.stim.visDiff==0 & uniBlk.tri.stim.audDiff<0);
else
    contra = uniBlk.tri.stim.visDiff>0 | (uniBlk.tri.stim.visDiff==0 & uniBlk.tri.stim.audDiff>0);
end
% subset = find(contra==1);
% subset = subset(randperm(length(subset)));
% contra = contra*0;
% contra(subset(1:400)) = 1;
% xlim([xDat(1) xDat(end)])
% disp(sum(contra))
for i = 1:2
    tDat.fracR(i,:) = cellfun(@(x) mean(uniBlk.tri.inactivation.data2Use(x & contra)), sampIdx);
    tDat.nTri(i,:) = cellfun(@(x) length(uniBlk.tri.inactivation.data2Use(x & contra)), sampIdx);
    tDat.nRT(i,:) = cellfun(@(x) sum(uniBlk.tri.inactivation.data2Use(x & contra)), sampIdx);
end

if plotType==1
    for i = 1
        CI = (1.96*sqrt((tDat.fracR(i,:).*(1-tDat.fracR(i,:)))./tDat.nTri(i,:)));
        pltM = tDat.fracR(i,:)-nrmDat.fracR;
        plotData = cat(3, pltM, pltM-CI, pltM+CI);
        opt.Marker = 'none';
        plt.rowsOfGrid(xDat', plotData, pltCol/i, opt);
        
        for j = 1:length(xDat)
            tbl = table([nrmDat.nRT;tDat.nRT(i,j)],[nrmDat.nTri-nrmDat.nRT;tDat.nTri(i,j)-tDat.nRT(i,j)], ...
                'VariableNames',{'RT','LT'},'RowNames',{'Cnt','Las'});
            [~,pVal(1,j)] = fishertest(tbl);
%             pVal(1,j) = prc.chiTest([nrmDat.nRT;tDat.nRT(i,j)], [nrmDat.nTri;tDat.nTri(i,j)]);
        end
        if any(pVal<pLim); plot(xDat(pVal<pLim), max(pltM+CI)+0.05, '.', 'color',  pltCol/i); end
    end
else
    CI = 1.96*sqrt((tDat.fracR(1,:).*(1-tDat.fracR(1,:)))./tDat.nTri(1,:));
    pltM = diff(tDat.fracR)*-1;
    plotData = cat(3, pltM, pltM-CI, pltM+CI);
    opt.Marker = 'none';
    plt.rowsOfGrid(xDat', plotData, pltCol, opt);
    
    pVal = runningmax<pLim
    
    for j = 1:length(xDat)
        pVal(1,j) = prc.chiTest([tDat.nRT(1,j);tDat.nRT(2,j)], [tDat.nTri(1,j);tDat.nTri(2,j)]);
    end
    if any(pVal<pLim); plot(xDat(pVal<pLim), max(ylim), '.', 'color',  pltCol/i); end
end
ylim(yRng)
plot(xlim, [0,0], '--k', 'linewidth', 1.5)
plot([0,0], ylim, '--k', 'linewidth', 1.5)
end
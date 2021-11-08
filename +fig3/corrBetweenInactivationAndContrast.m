function corrBetweenInactivationAndContrast(opt)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
uniBlks = spatialAnalysis('all', 'uniscan', 0, 1);
%%
if ~exist('opt', 'var'); opt = struct; end
if ~isfield(opt, 'prmType'); opt.prmType = 'rea'; end
if ~isfield(opt, 'pltType'); opt.pltType = 'mean'; end
if ~isfield(opt, 'revNorm'); opt.revNorm = 1; end
if ~isfield(opt, 'siteLoc'); opt.siteLoc = 'contra'; end
if ~isfield(opt, 'offset'); opt.offset = 1; end
%pre-assign performance and reaction structures with nans
nMice = length(uniBlks.blks);
clear reacN reacL
galvoIdx = {[0.6 2; 1.8, 2; 0.6, 3];[1.8 -4; 3,-4; 3,-3];[4.2,-2; 4.2,-3; 4.2,-4]};
[LME.reacT, LME.stimC, LME.mIdx, LME.lasOn] = deal(cell(length(galvoIdx), 2, 2));

for mouse = 1:nMice
    iBlk = prc.filtBlock(uniBlks.blks(mouse), uniBlks.blks(mouse).tri.inactivation.galvoPosition(:,2)~=4.5);
    iBlk = prc.filtBlock(iBlk, ~ismember(abs(iBlk.tri.inactivation.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | iBlk.tri.inactivation.laserType==0);
    iBlk = prc.filtBlock(iBlk, iBlk.tri.trialType.repeatNum==1 & iBlk.tri.trialType.validTrial & ~iBlk.tri.trialType.blank);
    
    if strcmpi(opt.prmType, 'rea')
        f2Use = 'reactionTimeComb';
        op2use  = @nanmedian;
        iBlk.tri.datGlobal = iBlk.tri.outcome.reactionTime;
        iBlk = prc.filtBlock(iBlk, ~isnan(iBlk.tri.outcome.reactionTime));
    else
        f2Use = 'fracTimeOutComb';
        op2use  = @mean;
        iBlk.tri.datGlobal = iBlk.tri.outcome.responseRecorded==0;
    end
    
    idx2Flip = iBlk.tri.inactivation.galvoPosition(:,1)<0 & iBlk.tri.inactivation.laserType==1;
    iBlk.tri.stim.audDiff(idx2Flip) = -1*iBlk.tri.stim.audDiff(idx2Flip);
    iBlk.tri.stim.visDiff(idx2Flip) = -1*iBlk.tri.stim.visDiff(idx2Flip);
    iBlk.tri.inactivation.galvoPosition(idx2Flip,1) = -1*iBlk.tri.inactivation.galvoPosition(idx2Flip,1);
    iBlk.tri.stim.conditionLabel(idx2Flip) = -1*iBlk.tri.stim.conditionLabel(idx2Flip);
    
    tBlk = iBlk;
    if opt.revNorm
        rIdx = (tBlk.tri.stim.visDiff>0 | (tBlk.tri.stim.visDiff==0 & tBlk.tri.stim.audDiff>0)) & tBlk.tri.inactivation.laserType==0;
        tBlk.tri.stim.audDiff(rIdx) = tBlk.tri.stim.audDiff(rIdx)*-1;
        tBlk.tri.stim.visDiff(rIdx) = tBlk.tri.stim.visDiff(rIdx)*-1;
        tBlk.tri.stim.conditionLabel(rIdx) = -1*tBlk.tri.stim.conditionLabel(rIdx);
    end
    
    rIdx = (tBlk.tri.stim.visDiff>0 | (tBlk.tri.stim.visDiff==0 & tBlk.tri.stim.audDiff>0));
    if strcmpi(opt.siteLoc, 'contra')
        tBlk = prc.filtBlock(tBlk, ~rIdx);
    else
        tBlk = prc.filtBlock(tBlk, rIdx);
    end
    tTypeIdx = {tBlk.tri.trialType.conflict; tBlk.tri.trialType.conflict};
    
    for site = 1:length(galvoIdx)
        gIdx = ismember(tBlk.tri.inactivation.galvoPosition,galvoIdx{site}, 'rows');
        for triT = 1:length(tTypeIdx)
            fBlk = tBlk;
            fBlk = prc.filtBlock(fBlk, (fBlk.tri.inactivation.laserType==0 | gIdx) & tTypeIdx{triT});
            normBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==0);
            lasBlk = prc.filtBlock(fBlk, fBlk.tri.inactivation.laserType==1);
            
            lasGrds = prc.getGridsFromBlock(lasBlk, 3);
            normGrds = prc.getGridsFromBlock(normBlk, 3);
            
            tkIdx = ~isnan(lasGrds.(f2Use));
            LME.reacT{site,triT,2} = [LME.reacT{site,triT,2}; lasGrds.(f2Use)(tkIdx)];
            LME.stimC{site,triT,2} = [LME.stimC{site,triT,2}; lasGrds.visValues(tkIdx)];
            LME.mIdx{site,triT,2} = [LME.mIdx{site,triT,2}; lasGrds.(f2Use)(tkIdx)*0+mouse];
            
            tkIdx = ~isnan(normGrds.(f2Use));
            LME.reacT{site,triT,1} = [LME.reacT{site,triT,1}; normGrds.(f2Use)(tkIdx)];
            LME.stimC{site,triT,1} = [LME.stimC{site,triT,1}; normGrds.visValues(tkIdx)];
            LME.mIdx{site,triT,1} = [LME.mIdx{site,triT,1}; normGrds.(f2Use)(tkIdx)*0+mouse];
            disp([site triT]);
        end
    end
end
%%
stimC = abs(LME.stimC{1,1,1});
deltaRT = LME.reacT{1,1,1}-LME.reacT{1,2,2};

LMEtlbs = table(deltaRT,stimC,LME.mIdx{1,1,1}, 'VariableNames',{'ReactTime','vCont','Mouse'});
LMEfits = fitlme(LMEtlbs, 'ReactTime~vCont+(1|Mouse)');
LMEVal = LMEfits.Coefficients.Estimate;


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
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 200 figWidth, figHeight]);

plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;

deltaRT = arrayfun(@(x) LME.reacT{1,1,1}(LME.mIdx{1,1,1}==x)-LME.reacT{1,2,2}(LME.mIdx{1,1,1}==x), 1:5, 'uni', 0);
stimC = arrayfun(@(x) abs(LME.stimC{1,1,1}(LME.mIdx{1,1,1}==x)), 1:5, 'uni', 0);
hold on
for i = 1:5
    plot(stimC{i}, deltaRT{i}, 'color', [0.5 0.5 0.5])
end
plot(stimC{1}, stimC{1}.*LMEVal(2)+LMEVal(1), 'k', 'linewidth', 2)
%%
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\corrBetweenInactivationAndContrast', '-pdf', '-painters');


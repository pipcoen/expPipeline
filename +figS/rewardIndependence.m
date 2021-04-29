function behaviourDetails(behBlks)
%% This function plots the data panels for figure one of the ms
if ~exist('behBlks', 'var') || isempty(behBlks); behBlks = spatialAnalysis('all', 'behavior', 1, 0); end

%%
nBlk = spatialAnalysis.getBlockType(behBlks.blks, 'norm');
nBlk = prc.filtBlock(nBlk, ~(nBlk.tri.stim.visContrast==0.06));
vValues = [-1*[0.8 0.4 0.2 0.1] 0 0.1 0.2 0.4 0.8];
aValues = [-60 0 60];
[vGrid, aGrid] = meshgrid(vValues,aValues);

numTrials = arrayfun(@(x,y) sum(nBlk.tri.stim.visDiff==x & nBlk.tri.stim.audDiff==y), vGrid, aGrid);
pRRewarded = double(vGrid > 0 & aGrid > 0 | vGrid > 0 & aGrid == 0 | aGrid > 0 & vGrid == 0);
pRRewarded(vGrid == 0 & aGrid == 0 | vGrid.*aGrid < 0) = 0.5;
probTrialPresented = numTrials./sum(numTrials(:));

pVAGivenR = probTrialPresented.*pRRewarded*2;
margResult = sum(pVAGivenR,2)*sum(pVAGivenR,1);

figure;
subplot(2,1,1);
imagesc(pVAGivenR); colormap gray;
axis equal;
box off;
set(gca,'TickLength',[0 0])
xticks(1:9);
set(gca, 'xticklabels', vValues)

subplot(2,1,2);
imagesc(margResult); colormap gray;
axis off
axis equal;
end
function movementPlots(behBlks)
%%
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1, 'raw'); end

%%
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
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);

%%

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
behBlks.viewRightLeftWheelSeparationOverTime;
axis square; 
xlim([0 0.3]);
box off;

axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
behBlks.viewConflictChoiceVsRT('kss');
axis square; 
xlim([0 0.3]);
box off;

axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);
behBlks.viewConflictChoiceVsRT('ksm');
axis square; 
xlim([0 0.3]);
box off;

axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg);
behBlks.viewConflictChoiceVsRT('fra');
axis square; 
xlim([0 0.3]);
ylim([0.3 0.7]);
box off;

export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\1_movementPlots', '-pdf', '-painters');
end
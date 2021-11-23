function timedInactivations
s = spatialAnalysis('all', 'variscan', 1, 1);

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
conT = 1;

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('v1', 'stim', 'vis',1,conT);

% axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('mos', 'stim', 'vis',1,conT);

axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('v1', 'stim', 'aud',1,conT);

% axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('mos', 'stim', 'aud',1,conT);
%%
% export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\timedInactivations', '-pdf', '-painters');
%%
% 
% 
% 
% axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla;
% s.viewInactivationTimedEffects('v1', 'move',1);
% 
% axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); cla;
% s.viewInactivationTimedEffects('mos', 'move',1);
% 
% axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); cla;
% s.viewInactivationTimedEffects('v1', 'move',2);
% 
% axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg); cla;
% s.viewInactivationTimedEffects('mos', 'move',2);

%%
end
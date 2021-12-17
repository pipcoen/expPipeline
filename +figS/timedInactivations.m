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
contra = 1;
yRng = [-0.1 0.4];
pType = 'res';

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('v1', 'vis', pType,contra, yRng);
s.viewInactivationTimedEffects('v1', 'aud', pType,contra, yRng);
title('VIS Inactivation')

axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla;
s.viewInactivationTimedEffects('mos', 'vis', pType,contra, yRng);
s.viewInactivationTimedEffects('mos', 'aud', pType,contra, yRng);
title('MOs Inactivation')

axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); cla;
yRng = [-0.2 0.5];
s.viewInactivationTimedEffects('s1', 'vis', pType,contra, yRng);
s.viewInactivationTimedEffects('s1', 'aud', pType,contra, yRng);
title('S1 Inactivation')

axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); cla;
yRng = [-0.2 0.4];
s.viewInactivationTimedEffects('out', 'vis', pType,contra, yRng);
s.viewInactivationTimedEffects('out', 'aud', pType,contra, yRng);
title('Outside Brain')

plt.suplabel('Change in Rightward Trials with Contra Inactivation', 't')

% export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\timedInactivations_Choices_Ipsi', '-pdf', '-painters');
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
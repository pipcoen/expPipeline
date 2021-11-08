function mouseWith5AudLevels
%% This function plots the data panels for figure one of the ms
s = spatialAnalysis('all', 'aud5', 1, 1);
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
%%
axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);cla
s.viewGLMFits('simpLogSplitVEmpA', [],[], 1)
axis square;

%%
axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg);cla
s.viewGLMFits('simpLogSplitVEmpA', [],'log', 1)
axis square;
ylim([-3 3]);

%%
axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);cla
s.viewGLMFits('simpEmp', [],[], 1)
axis square;

%%
axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg);cla
s.viewGLMFits('simpEmp', [],'log', 1)
axis square;
ylim([-3 3]);

%%
s.viewGLMFits('simpLogSplitVEmpA', 5,'none', 1)
axis square;
LL{3,1} = s.glmFit{1}.logLik;

%%
s.viewGLMFits('fullEmp', 5,'none', 1)
axis square;
LL{4,1} = s.glmFit{1}.logLik; 
%%
s.viewGLMFits('visOnly', 5,'none', 1)
axis square;
LL{1,1} = s.glmFit{1}.logLik; 

s.viewGLMFits('audOnly', 5,'none', 1)
axis square;
LL{2,1} = s.glmFit{1}.logLik; 

%%
s.viewGLMFits('biasOnly', 5,'none', 1)
offset = mean(s.glmFit{1}.logLik); 

%%
axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);cla
plt.jitter(cellfun(@(x) offset-x, LL, 'uni', 0));
set(gca, 'position', get(gca, 'position').*[1 1 0.5 1])
xlim([0.8 4.1])
ylim([0 0.5])
set(gca, 'xTickLabel', {'visOnly', 'audOnly', 'Add', 'Full'})

[~,~,stats] = anova1(cell2mat(LL'));
%%
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\mouseWith5AudLevels', '-pdf', '-painters');
end
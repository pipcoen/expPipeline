function plotOpt = compareModels(fileNames, minTrials)
% if ~exist('fileNames', 'var'); fileNames = {'CrossVal_AudOnly';'CrossVal_VisOnly'; 'CrossVal_AudDom'; 'CrossVal_SimpLog'}; end
% if ~exist('fileNames', 'var'); fileNames = {'CrossVal_ReducedLogCN'}; end
% if ~exist('fileNames', 'var'); fileNames = {'simpLog_Cross5';'simpLogSplitA_Cross5';'simpLogSplitV_Cross5'; 'simpLogSplitVSplitA_Cross5'}; end
if ~exist('fileNames', 'var'); fileNames = {'audOnly_Cross5';'visOnly_Cross5';'simpLogSplitVSplitA_Cross5'}; end
if ~exist('minTrials', 'var'); minTrials = 1; end
plotOpt.yData = cell(length(fileNames),1);
if any(contains(fileNames, 'Nest')); load FullEmpNestMaxPerformance s;
else, load FullEmpMaxPerformance.mat s;
end
included = cellfun(@(x) sum(~isnan(x.blockData.tri.outcome.responseCalc)), s.glmFit)>minTrials;
maxPerformance(:,1) = arrayfun(@(x) cell2mat(unique(x.exp.subject)'), s.blks(included), 'uni', 0);
maxPerformance(:,2) = cellfun(@(x) x.logLik, s.glmFit(included), 'uni', 0)';

if any(contains(fileNames, 'Nest')); load BiasOnly4TONestPerformance s;
else, load BiasOnlyPerformance.mat s;
end
minPerformance(:,1) = arrayfun(@(x) cell2mat(unique(x.exp.subject)'), s.blks(included), 'uni', 0);
minPerformance(:,2) = cellfun(@(x) x.logLik, s.glmFit(included), 'uni', 0)';

plotOpt.xTickLabels = ['Bias'; fileNames; 'Full'];
plotOpt.yData = cell(length(plotOpt.xTickLabels),1);
plotOpt.yData{1} = cell2mat(minPerformance(:,2))*-1;
plotOpt.yData{end} = cell2mat(maxPerformance(:,2))*-1;
plotOpt.subjects = minPerformance(:,1);

for i = 1:length(fileNames)
    load(fileNames{i});
    logLik = cell2mat(cellfun(@(x) mean(x.logLik), s.glmFit(included), 'uni', 0));
    plotOpt.yData{i+1} = logLik*-1;
end
expType = arrayfun(@(x) cell2mat(unique(x.exp.expType)'), s.blks, 'uni', 0);
plotOpt.faceColors = {repmat([0 0 0], length(expType), 1)};
plotOpt.faceColors{1}(strcmp(expType, 'inactivation'),:) = repmat([0 1 1], sum(strcmp(expType, 'inactivation')),1);
plotOpt.faceColors = repmat(plotOpt.faceColors, length(plotOpt.yData),1);
figure;
plt.jitter(plotOpt.yData, plotOpt); grid('on');
end

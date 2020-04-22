function compareModels(fileNames, minTrials)
% if ~exist('fileNames', 'var'); fileNames = {'CrossVal_AudOnly';'CrossVal_VisOnly'; 'CrossVal_AudDom'; 'CrossVal_SimpLog'}; end
% if ~exist('fileNames', 'var'); fileNames = {'CrossVal_ReducedLogCN'}; end
if ~exist('fileNames', 'var'); fileNames = {'simpLog_Cross5';'simpLogSplitA_Cross5';'simpLogSplitV_Cross5'; 'simpLogSplitVSplitA_Cross5'}; end
if ~exist('minTrials', 'var'); minTrials = 1; end
plotOpt.yData = cell(length(fileNames),1);
if any(contains(fileNames, 'Nest')); load FullEmpNestMaxPerformance s;
else, load FullEmpMaxPerformance.mat s; 
end
included = cellfun(@(x) length(x.blockData.tri.outcome.responseMade), s.glmFit)>minTrials;
maxPerformance(:,1) = arrayfun(@(x) unique(x.exp.subject), s.blks(included));
maxPerformance(:,2) = cellfun(@(x) x.logLik, s.glmFit(included), 'uni', 0)';

if any(contains(fileNames, 'Nest')); load BiasOnly4TONestPerformance s;
else, load BiasOnlyPerformance.mat s; 
end
minPerformance(:,1) = arrayfun(@(x) unique(x.exp.subject), s.blks(included));
minPerformance(:,2) = cellfun(@(x) x.logLik, s.glmFit(included), 'uni', 0)';
plotOpt.xTickLabels = fileNames;
for i = 1:length(fileNames)
    load(fileNames{i});  
    logLik = cell2mat(cellfun(@(x) mean(x.logLik), s.glmFit(included), 'uni', 0));
    plotOpt.yData{i,1} = 100*(cell2mat(minPerformance(:,2))-logLik)./(cell2mat(minPerformance(:,2))-cell2mat(maxPerformance(:,2)));
end
expType = arrayfun(@(x) unique(x.exp.expType), s.blks);
plotOpt.faceColors = {repmat([0 0 0], length(expType), 1)};
plotOpt.faceColors{1}(strcmp(expType, 'inactivation'),:) = repmat([0 1 1], sum(strcmp(expType, 'inactivation')),1);
plotOpt.faceColors = repmat(plotOpt.faceColors, length(plotOpt.yData),1);
% plotOpt.pairs2test = 'all';
figure;
plt.jitter(plotOpt.yData, plotOpt); grid('on');
ylim([90 100])
end

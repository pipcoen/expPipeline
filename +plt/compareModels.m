function compareModels(fileNames, minTrials)
if ~exist('minTrials', 'var'); minTrials = 1; end
plotOpt.yData = cell(length(fileNames),1);
if any(contains(fileNames, 'Nest')); load FullEmpNestMaxPerformance s;
else, load FullEmpMaxPerformance.mat s; 
end
included = cellfun(@(x) length(x.blockData.responseMade), s.glmFit)>minTrials;
maxPerformance(:,1) = s.subjects(included);
maxPerformance(:,2) = cellfun(@(x) x.logLik, s.glmFit(included), 'uni', 0)';

if any(contains(fileNames, 'Nest')); load BiasOnly4TONestPerformance s;
else, load FullEmpMaxPerformance.mat s; 
end
minPerformance(:,1) = s.subjects(included);
minPerformance(:,2) = cellfun(@(x) x.logLik, s.glmFit(included), 'uni', 0)';

plotOpt.xTickLabels = fileNames;
for i = 1:length(fileNames)
    load(fileNames{i});  
    logLik = cell2mat(cellfun(@(x) mean(x.logLik), s.glmFit(included), 'uni', 0))';
    plotOpt.yData{i,1} = 100*(cell2mat(minPerformance(:,2))-logLik)./(cell2mat(minPerformance(:,2))-cell2mat(maxPerformance(:,2)));
end
figure;
plt.jitter(plotOpt.yData, plotOpt); grid('on');
end

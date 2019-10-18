function viewLearningRate(subjects)
figure;
axesOpt.totalNumOfAxes = 3;
axesOpt.btlrMargins =  [80 100 80 40];
axesOpt.gapBetweenAxes = [100 60];
numSessions = 40;
multiTest = 20;
tempObj = spatialAnalysis(unique(subjects), ['first' num2str(numSessions)], -1, 'none', 'multiSpaceWorld');
%%
allPerformance = arrayfun(@(x) cell2mat(x.exp.performanceAVM),tempObj.blks, 'uni', 0);
indiColor = [[0.5 0 0.5]; [0.9294 0.6902 0.1294]/2; [0.5 0.5 0.5]];
meanColor = [[1 0 1]; [0.9294 0.6902 0.1294]; [0 0 0]];
tempObj.hand.axes = plt.getAxes(axesOpt,1); hold on
firstMulti = cell2mat(cellfun(@(x) x(1:multiTest,3)',allPerformance, 'uni', 0));
% plt.stdPatch(1:multiTest, nanmean(firstMulti), nanstd(firstMulti), indiColor(3,:), 0.3, 1)
plot(1:multiTest,firstMulti, 'color', indiColor(3,:));
plot(nanmean(firstMulti), 'color', meanColor(3,:), 'LineWidth', 3);
box off;
ylim([40 100]); xlim([1 multiTest])
plot([1 multiTest],[50, 50], '--k')
title(sprintf('First sessions from %d mice', length(subjects)));

tempObj.hand.axes = plt.getAxes(axesOpt,2); hold on
for i = 1:2
    uniIdx = num2cell(cellfun(@(x) find(~isnan(x(:,i)),1), allPerformance));
    first10Uni = cell2mat(cellfun(@(x,y) x(y:y+9,i)',allPerformance, uniIdx, 'uni', 0));
    plt.stdPatch(1:10, nanmean(first10Uni), nanstd(first10Uni), indiColor(i,:), 0.3, 1)
    plot(nanmean(first10Uni), 'color', meanColor(i,:), 'LineWidth', 3);
end
box off;
ylim([40 100]); xlim([1 multiTest])
plot([1 multiTest],[50 50], '--k')
title(sprintf('First unisensory sessions from %d mice', length(tempobj.blks.exp.subject)));
% 
% %%
% firstDay = cell2mat(arrayfun(@(x) datenum(x.params(1).expDate),tempObj.blks, 'uni', 0));
% tempObj = spatialAnalysis(obj.blks.exp.subject(~contains(obj.blks.exp.subject, 'DJ')), 'last');
% lastDay = cell2mat(cellfun(@(x) datenum(x.params(1).expDate),tempObj.blks, 'uni', 0));
% oldMice = (lastDay-firstDay)>200;
% dateRanges = arrayfun(@(x) {'rng', datestr(x+180, 'yyyy-mm-dd') datestr(x+180+numSessions, 'yyyy-mm-dd')}, firstDay, 'uni', 0);
% tempObj = spatialAnalysis(tempobj.blks.exp.subject(oldMice), dateRanges(oldMice));
% 
% %%
% tempObj.hand.axes = plt.getAxes(axesOpt,3); hold on
% allPerformance = cellfun(@(x) [[x.params.audPerformance]' [x.params.visPerformance]' [x.params.mulPerformance]'],tempObj.blks, 'uni', 0);
% for i = 1:3
%     uniIdx = num2cell(cellfun(@(x) find(~isnan(x(:,i)),1), allPerformance));
%     first10Uni = cell2mat(cellfun(@(x,y) x(y:y+9,i)',allPerformance, uniIdx, 'uni', 0));
%     plt.stdPatch(1:10, nanmean(first10Uni), nanstd(first10Uni), indiColor(i,:), 0.3, 1)
%     plot(nanmean(first10Uni), 'color', meanColor(i,:), 'LineWidth', 3);
% end
% box off;
% ylim([40 100]); xlim([1 10])
% plot([1 10],[50, 50], '--k')
% title(sprintf('10 Sessions after 6 months from %d mice', length(tempobj.blks.exp.subject)));
end
function viewLearningRate(obj)
figure;
axesOpt.totalNumOfAxes = 3;
axesOpt.btlrMargins =  [80 100 80 40];
axesOpt.gapBetweenAxes = [100 60];
allPerformance = cell(length(obj.subjects),3);
for i  = 1:length(obj.subjects)
    tempObj = copy(obj);
    tempObj = tempObj.changeMouse(tempObj.subjects(i), {'first200'}, 0, 'prm');
    tempDat = [tempObj.params{1}.audPerformance tempObj.params{1}.visPerformance tempObj.params{1}.mulPerformance];
    tempDat(tempObj.params{1}.validResponses<10,:) = [];
    allPerformance(i,:) = arrayfun(@(x) tempDat(~isnan(tempDat(:,x)),x),1:3,'uni', 0);
end

indiColor = [[0 0.5 0.5]; [0.5 0 0.5]; [0.5 0.5 0.5]];
meanColor = 'cmk';
obj.axesHandles = plt.getAxes(axesOpt,1); hold on
cellfun(@(x) plot(x(1:20), 'color', indiColor(3,:)),allPerformance(:,3));
plot(mean(cell2mat(cellfun(@(x) x(1:20), allPerformance(:,3)', 'uni', 0)),2), meanColor(3), 'linewidth', 3)
box off;
ylim([40 100]); xlim([1 20])
plot([1 20],[50, 50], '--k')

obj.axesHandles = plt.getAxes(axesOpt,2); hold on
ylim([40 100]); xlim([1 20])
for i = 1:2
    cellfun(@(x) plot(x(1:20), 'color', indiColor(i,:)),allPerformance(:,i));
    plot(mean(cell2mat(cellfun(@(x) x(1:20), allPerformance(:,i)', 'uni', 0)),2), meanColor(i), 'linewidth', 3)
end
box off;
plot([1 20],[50 50], '--k')

obj.axesHandles = plt.getAxes(axesOpt, 3); hold on;
allPerformance = allPerformance(cellfun(@length, allPerformance(:,3))>100,:);
for i = 1:3
    cellfun(@(x) plot(x(end-19:end), 'color', indiColor(i,:)),allPerformance(:,i));
    plot(mean(cell2mat(cellfun(@(x) x(end-19:end), allPerformance(:,i)', 'uni', 0)),2), meanColor(i), 'linewidth', 3)
end
box off;
ylim([40 100]); xlim([1 20])
plot([1 20],[50, 50], '--k')
end
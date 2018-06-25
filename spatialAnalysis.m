function viewLearningRate(obj)
figure;
allPerformance = nan*ones(500,3,length(obj.subjects));
axesOpt.totalNumOfAxes = 3;
axesOpt.btlrMargins =  [80 100 80 40];
axesOpt.gapBetweenAxes = [100 60];
axesOpt.idx = 1;
for i  = 1:length(obj.subjects)
    tempObj = copy(obj);
    tempObj = tempObj.changeMouse(tempObj.subjects(i), {'first200'}, 0, 'prm');
    tempDat = [tempObj.params{1}.audPerformance tempObj.params{1}.visPerformance tempObj.params{1}.mulPerformance];
    tempDat(tempObj.params{1}.validResponses<10,:) = [];
    srtIdx = find(~isnan(tempObj.params{1}.mulPerformance),1);
    allPerformance(1:size(tempDat,1)-srtIdx+1, :, i) = tempDat;
end
obj.axesHandles = plt.getAxes(axesOpt);
plot(1:20, squeeze(allPerformance(1:20,3,:))', 'color', [0.5 0.5 0.5]); hold on;
plot(1:20, mean(squeeze(allPerformance(1:20,3,:)),2), 'k', 'linewidth', 3)
box off;
ylim([40 100]); xlim([1 20])
axesOpt.idx = 2;
obj.axesHandles = plt.getAxes(axesOpt);
for i  = 1:length(obj.subjects)
    tDat = allPerformance(~any(isnan(allPerformance(:,:,i)),2),:,i);
    initialAudVis(1:25,:,i) = tDat(1:25,:);
    for j = 1:3
        latestAudVis(1:10,j,i) = smooth(tDat(end-9:end,j),1);
    end
end
ylim([40 100]); xlim([1 25])

indiColor = [[0 0.5 0.5]; [0.5 0 0.5]; [0.5 0.5 0.5]];
mCol = 'cmk';
for i = 1:3
    plot(1:size(initialAudVis,1), squeeze(initialAudVis(:,i,:))', 'color', indiColor(i,:)); hold on;
    plot(1:size(initialAudVis,1), mean(squeeze(initialAudVis(:,i,:)),2), mCol(i), 'linewidth', 3)
    box off;
end
ylim([40 100]); xlim([1 25])
axesOpt.idx = 3;
obj.axesHandles = plt.getAxes(axesOpt);
for i = 1:3
    plot(1:size(latestAudVis,1), squeeze(latestAudVis(:,i,3:5))', 'color', indiColor(i,:)); hold on;
    plot(1:size(latestAudVis,1), mean(squeeze(latestAudVis(:,i,3:5)),2), mCol(i), 'linewidth', 3)
    box off;
end
ylim([40 100]); xlim([1 50])
end
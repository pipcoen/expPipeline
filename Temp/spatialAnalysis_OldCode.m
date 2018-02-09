function viewLearningRate(obj)
            figure;
            allPerformance = nan*ones(500,3,length(obj.subjects));
            for i  = 1:length(obj.subjects)
                tempObj = obj.changeMouse(obj.subjects(i), {'all'}, 'prm');
                numSessions = length(tempObj.params{1}.mulPerformance);
                srtIdx = find(~isnan(tempObj.params{1}.mulPerformance),1);
                allPerformance(1:numSessions-srtIdx+1, 1, i) = tempObj.params{1}.audPerformance(srtIdx:end);
                allPerformance(1:numSessions-srtIdx+1, 2, i) = tempObj.params{1}.visPerformance(srtIdx:end);
                allPerformance(1:numSessions-srtIdx+1, 3, i) = tempObj.params{1}.mulPerformance(srtIdx:end);
            end
            plt.getAxes(1, 3, [], [], [80 100 80 40], [100 60]);
            plot(1:25, squeeze(allPerformance(1:25,3,:))', 'color', [0.5 0.5 0.5]); hold on;
            plot(1:25, mean(squeeze(allPerformance(1:25,3,:)),2), 'k', 'linewidth', 3)
            box off;
            ylim([40 100]); xlim([1 25])
            plt.getAxes(2, 3, [], [], [80 100 80 40], [100 60]);
            for i  = 1:length(obj.subjects)
                tDat = allPerformance(~any(isnan(allPerformance(:,:,i)),2),:,i);
                initialAudVis(1:25,:,i) = tDat(1:25,:);
                for j = 1:3
                    latestAudVis(1:50,j,i) = smooth(tDat(end-49:end,j),1);
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
            
            plt.getAxes(3, 3, [], [], [80 100 80 40], [100 60]);
            for i = 1:3
                plot(1:size(latestAudVis,1), squeeze(latestAudVis(:,i,3:5))', 'color', indiColor(i,:)); hold on;
                plot(1:size(latestAudVis,1), mean(squeeze(latestAudVis(:,i,3:5)),2), mCol(i), 'linewidth', 3)
                box off;
            end
            ylim([40 100]); xlim([1 50])
        end
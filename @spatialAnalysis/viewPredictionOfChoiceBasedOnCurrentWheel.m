function viewPredictionOfChoiceBasedOnCurrentWheel(obj)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices
%	'rea'--------------------timeToFirstMove vs fration of rightward choices

% xDat = 0:50:500;
xDat = 0:0.005:0.2;
tStep = xDat(2)-xDat(1);
moveDatAVC = cell(3,2);
colorBuff = floor(length(xDat)/20)*2+1;
figure;
hold on;
for i  = 1:length(obj.blks)
    bigBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    bigBlk = prc.filtBlock(bigBlk, ~isnan(bigBlk.tri.outcome.timeToFirstMove) & bigBlk.tri.outcome.timeToFirstMove<0.5);
    bigBlk.tri.right = (bigBlk.tri.stim.visInitialAzimuth>0 & bigBlk.tri.trialType.visual) | bigBlk.tri.stim.audInitialAzimuth > 0;
    bigBlk.tri.left = (bigBlk.tri.stim.visInitialAzimuth<0 & bigBlk.tri.trialType.visual) | bigBlk.tri.stim.audInitialAzimuth < 0;
    contast2Use = max(abs(bigBlk.tri.stim.visContrast(bigBlk.tri.trialType.visual)));
    bigBlk = prc.filtBlock(bigBlk, (abs(bigBlk.tri.stim.visContrast) == contast2Use) | bigBlk.tri.trialType.auditory);
    rawWheelTV = bigBlk.tri.raw.wheelTimeValue;

    %Align to first move
    alignTimes = repmat({xDat}, bigBlk.tot.trials,1);
%     alignTimes = arrayfun(@(x) x+xDat, bigBlk.tri.outcome.timeToFirstMove, 'uni', 0);
%     
    bigBlk.tri.wheelPos = cellfun(@(x,y) interp1(x(:,1), x(:,2), [y y(end)+tStep], 'nearest', 'extrap')*-1, rawWheelTV, alignTimes, 'uni', 0);
    bigBlk.tri.wheelVel = cellfun(@diff, bigBlk.tri.wheelPos, 'uni', 0);
    bigBlk.tri.wheelPos = cellfun(@(x) x(1:end-1), bigBlk.tri.wheelPos, 'uni', 0);
    
%     filtIdx = bigBlk.tri.trialType.visual;
%     bigBlk = prc.filtBlock(bigBlk, filtIdx);
%     alignTimes = alignTimes(filtIdx);
        
    colorMap = sqrt(plt.redBlueMap(length(xDat)*2+colorBuff));
    colors2Use{1,1} = flipud(colorMap(1:length(xDat),:));
    colors2Use{2,1} = colorMap(end-length(xDat)+1:end,:);
    [meanPos, meanVel, sqrDis] = deal(cell(2,1));
    
    threshMoveTime = bigBlk.tri.outcome.threshMoveTime;
    allWheelPosRaw = cell2mat(arrayfun(@(x,y,z) [x{1}(z{1}<=y) nan*x{1}(z{1}>y)], bigBlk.tri.wheelPos, threshMoveTime, alignTimes, 'uni', 0));
    allWheelVelRaw = cell2mat(arrayfun(@(x,y,z) [x{1}(z{1}<=y) nan*x{1}(z{1}>y)], bigBlk.tri.wheelVel, threshMoveTime, alignTimes,  'uni', 0));
%     
%     allWheelPosRaw = cell2mat(bigBlk.tri.wheelPos);
%     allWheelVelRaw = cell2mat(bigBlk.tri.wheelVel);
    allWheelPos = normalize(allWheelPosRaw);
    allWheelVel = normalize(allWheelVelRaw);
    for j = 1:2
        wheelPos = allWheelPos(bigBlk.tri.outcome.responseMade==j, :);
        wheelVel = allWheelVel(bigBlk.tri.outcome.responseMade==j, :);
        meanPos{j} = nanmean(wheelPos);
        meanVel{j} = nanmean(wheelVel);
        
        scatter(meanPos{j}, meanVel{j}, [], colors2Use{j}, 'filled')
        sqrDis{j} = sqrt(bsxfun(@minus , meanPos{j}, allWheelPos).^2 + bsxfun(@minus , meanVel{j}, allWheelVel).^2);
        sqrDis{j}(allWheelPosRaw == 0 | allWheelVelRaw == 0) = nan;
    end
    bigBlk.tri.nearestMove = ((sqrDis{1}-sqrDis{2})>0) - ((sqrDis{1}-sqrDis{2})<0);
    
    filtIdx = {bigBlk.tri.trialType.auditory; bigBlk.tri.trialType.visual; bigBlk.tri.trialType.conflict};
    for j = 1:3
        subBlk = prc.filtBlock(bigBlk, filtIdx{j});
        nearestMove = subBlk.tri.nearestMove;
        if j == 2
            moveDatAVC{j,1}(i,:) = nanmean(nearestMove(subBlk.tri.right,:));
            moveDatAVC{j,2}(i,:) = nanmean(nearestMove(subBlk.tri.left,:));
        else,
            moveDatAVC{j,2}(i,:) = nanmean(nearestMove(subBlk.tri.right,:));
            moveDatAVC{j,1}(i,:) = nanmean(nearestMove(subBlk.tri.left,:));
        end
        moveDatAVC{j,3}(i,:) = moveDatAVC{j,1}(i,:) - moveDatAVC{j,2}(i,:);
    end
end
xlim([-1 1]);
ylim([-1 1]);
axis square
plot([-1 1], [0 0], '--k', 'linewidth', 3)
plot([0 0], [-1 1], '--k', 'linewidth', 3)
xlabel('Nomalized Displacement');
ylabel('Nomalized Velocity');
%%
% titles = {'Left Aud'; 'Right Aud'; 'Difference'};
% xLimEnd = find(xDat >= 0.2,1);
% moveDatAVCTuncated = cellfun(@(x) x(:,1:xLimEnd), moveDatAVC, 'uni', 0);
% xlabel('Time from stimulus onset (s)')
% ylabel('Similarity to rightward movement')
% for i = 1:3
%     meanDataAVC = cell2mat(cellfun(@(x) mean(x), moveDatAVCTuncated(:,i), 'uni', 0));
%     stdDataAVC = cell2mat(cellfun(@(x) std(x), moveDatAVCTuncated(:,i), 'uni', 0));
%     if i == 3; meanDataAVC = meanDataAVC/2; stdDataAVC = stdDataAVC/2; end
%     lowBoundAVC = meanDataAVC-stdDataAVC./sqrt(length(obj.blks));
%     upBoundAVC = meanDataAVC+stdDataAVC./sqrt(length(obj.blks));
%     
%     
%     figure;
%     plotData = cat(3, meanDataAVC, lowBoundAVC, upBoundAVC);
%     plotOpt.Marker = 'none';
%     colors = [1 0 1; 0 1 1; 0 0 0];
%     plt.rowsOfGrid(xDat(1:xLimEnd), plotData(:,:,:), colors, plotOpt);
%     plot([0 xDat(xLimEnd)], [0 0], '--k');
%     ylim([-0.4 0.4]);
%     title(titles{i});
% end
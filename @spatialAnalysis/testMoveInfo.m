function testMoveInfo(obj, testTag)
if ~exist('testTag', 'var'); testTag = ''; end
if ~strcmpi(testTag, 'svm')
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
    
for i = 1:6
    axesHandle = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
    blk = obj.blks(randperm(length(obj.blks),1));
    blk = spatialAnalysis.getBlockType(blk,'norm',1);
    
    timeToFirstMove = blk.tri.outcome.timeToFirstMove;
    reactionTime = blk.tri.outcome.reactionTime;
    timeToResponseThresh = blk.tri.outcome.timeToResponseThresh;
    
    idx = find(reactionTime~=timeToFirstMove & ~isnan(reactionTime.*timeToFirstMove));
    idx = idx(randperm(length(idx),1));
    
    timeToFirstMove = timeToFirstMove(idx);
    reactionTime = reactionTime(idx);
    choiceTime = timeToResponseThresh(idx);
    xDat = -0.05:0.001:choiceTime+0.005;
    rawWheelTV = blk.tri.raw.wheelTimeValue{idx};
    wheelPos = interp1(rawWheelTV(:,1), rawWheelTV(:,2), xDat, 'nearest', 'extrap');
    
    cla
    hold on;
    plot(xDat*1000, wheelPos, 'k');
    moveIdx = knnsearch(xDat', [timeToFirstMove; reactionTime; choiceTime]);
    plot(xDat(moveIdx(1))*1000, wheelPos(moveIdx(1)), '*r')
    plot(xDat(moveIdx(2))*1000, wheelPos(moveIdx(2)), '*b')
    plot(xDat(moveIdx(3))*1000, wheelPos(moveIdx(3)), '*m')
    plot([0 0], ylim, '--k');
end
end

%%
figure;
for i  = 1:length(obj.blks)
    numSamp = 1000;
    blk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    
    if strcmpi(testTag(1:3), 'svm')
        selectIdx = zeros(blk.tot.trials,1);
        
        lTri = find(cellfun(@(x) x(1,2), blk.tri.outcome.timeDirAllMoveOnsets) == 1);
        rTri = find(cellfun(@(x) x(1,2), blk.tri.outcome.timeDirAllMoveOnsets) == 2);
        
        selectIdx(lTri(randperm(length(lTri), round(numSamp/2)))) = 1;
        selectIdx(rTri(randperm(length(rTri), round(numSamp/2)))) = 1;
        blk = prc.filtBlock(blk, selectIdx);
        reactDir = blk.tri.outcome.responseCalc;
        reactTime = blk.tri.outcome.timeToFirstMove; 
%         reactTime = blk.tri.outcome.timeDirFirstMove(:,1); 
    end
    wheelSamples = -0.5:0.001:0.3;
    xIdx = arrayfun(@(x) x+(wheelSamples), reactTime, 'uni', 0);
    rawWheelTV = blk.tri.raw.wheelTimeValue;
    wheelPos = cell2mat(cellfun(@(x,y) interp1(x(:,1), x(:,2), y, 'linear', 'extrap'), rawWheelTV, xIdx, 'uni', 0));
    wheelVel = [wheelPos(:,1:10)*nan wheelPos(:,11:end)-wheelPos(:,1:end-10)];
    
    xAxis = [-0.15:0.01:0.1];    
    pnts2Take = ismember(round(wheelSamples*1e3), round(xAxis*1e3));
    if sum(pnts2Take)~=length(xAxis); error('SamplingError'); end 
    wheelPos = num2cell(wheelPos(:, pnts2Take),1);   
    wheelVel = num2cell(wheelVel(:, pnts2Take),1);
    
    svmMod.subject{i} = blk.exp.subject{1};
    svmMod.reactDir{i} = reactDir;
    svmMod.models(i,:) = cellfun(@(x) fitcsvm(x,reactDir, 'CrossVal', 'on', 'KFold', 4), wheelVel, 'uni', 0);
    svmMod.modPerf(i,:) = cell2mat(cellfun(@(x) mean(x.Y==x.kfoldPredict), svmMod.models(i,:), 'uni', 0));
    
    hold on;
    cCol = 'b';
    plot(xAxis,svmMod.modPerf(i,:), cCol);
    drawnow;
    disp(['Done: ' blk.exp.subject{1}]);
end
plot(xAxis,mean(svmMod.modPerf), cCol, 'linewidth',4);
% save svmData svmMod
end






%%


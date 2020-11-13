function testMoveInfo(obj, testTag)
if ~exist('testTag', 'var'); testTag = 'svm'; end
if ~strcmpi(testTag, 'svm')
for i = 1:200
    blk = obj.blks(randperm(length(obj.blks),1));
    blk = spatialAnalysis.getBlockType(blk,'norm',1);
    blk = prc.filtBlock(blk, blk.tri.outcome.timeDirChoiceInit<0.5);
    
    timeToFirstMove = blk.tri.outcome.timeDirFirstMove(:,1);
    timeToChoiceInit = blk.tri.outcome.timeDirChoiceInit(:,1);
    timeToChoiceCross = blk.tri.outcome.timeDirChoiceCross(:,1);
    
    idx = find(timeToChoiceInit==timeToFirstMove & ~isnan(timeToChoiceInit.*timeToFirstMove));
    idx = idx(randperm(length(idx),1));
    
    firstTime = timeToFirstMove(idx);
    choiceInitTime = timeToChoiceInit(idx);
    choiceTime = timeToChoiceCross(idx);
    xDat = -0.05:0.001:choiceTime+0.005;
    rawWheelTV = blk.tri.raw.wheelTimeValue{idx};
    wheelPos = interp1(rawWheelTV(:,1), rawWheelTV(:,2), xDat, 'linear', 'extrap');
    
    cla
    hold on;
    plot(xDat, wheelPos, 'k');
    moveIdx = knnsearch(xDat', [firstTime; choiceInitTime; choiceTime]);
    plot(xDat(moveIdx(1)), wheelPos(moveIdx(1)), '*r')
    plot(xDat(moveIdx(2)), wheelPos(moveIdx(2)), '*b')
    plot(xDat(moveIdx(3)), wheelPos(moveIdx(3)), '*m')
end
end



for i  = 1:length(obj.blks)
    numSamp = 1000;
    blk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    blk = prc.filtBlock(blk, blk.tri.outcome.timeDirChoiceInit<0.5);
    
    if strcmpi(testTag(1:3), 'svm')
        valTri = ~isnan(blk.tri.outcome.timeDirFirstMove(:,1)) & ~isnan(blk.tri.outcome.timeDirChoiceCross(:,1));
        selectIdx = zeros(blk.tot.trials,1);
        
        lTri = find(blk.tri.outcome.timeDirChoiceCross(:,2) == 1 & valTri);
        rTri = find(blk.tri.outcome.timeDirChoiceCross(:,2) == 2 & valTri);
        
        selectIdx(lTri(randperm(length(lTri), round(numSamp/2)))) = 1;
        selectIdx(rTri(randperm(length(rTri), round(numSamp/2)))) = 1;
        blk = prc.filtBlock(blk, selectIdx);
        reactDir = blk.tri.outcome.timeDirChoiceCross(:,2);
        reactTime = blk.tri.outcome.timeDirChoiceInit(:,1); 
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
    svmMod.models(i,:) = cellfun(@(x) fitcsvm(x,reactDir, 'CrossVal', 'on', 'KFold', 4), wheelPos, 'uni', 0);
    svmMod.modPerf(i,:) = cell2mat(cellfun(@(x) mean(x.Y==x.kfoldPredict), svmMod.models(i,:), 'uni', 0));
    
    hold on;
    plot(xAxis,svmMod.modPerf(i,:), 'b');
    drawnow;
    disp(['Done: ' blk.exp.subject{1}]);
end

save svmData svmMod
end






%%


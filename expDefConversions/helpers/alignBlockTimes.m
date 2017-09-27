function blkTConv = alignBlockTimes(block, timeline)
diodeChannel = timeline.hw.inputs(strcmp({timeline.hw.inputs.name},'photoDiode')).arrayColumn;
diodeSamples = timeline.rawDAQData(:,diodeChannel);
timeSamples = timeline.rawDAQTimestamps';
[~, diodeThreshold] = kmeans(diodeSamples, 2); diodeThreshold = mean(diodeThreshold);
diodeFlipPnts = abs(diff(diodeSamples > diodeThreshold)) > 0; % look for the flips
diodeFlipPnts = mean([timeSamples([false; diodeFlipPnts]) timeSamples([diodeFlipPnts; false])], 2);
blockFlipPnts = block.stimWindowUpdateTimes;

if length(blockFlipPnts)-length(diodeFlipPnts) > 10
    warning('Missed some flips on Photodiode. Will try to compensate');
    if mean(diodeFlipPnts(1:10) - blockFlipPnts(1:10)) < 0
        error('Could not compensate');
    end
    [covarValues, covarLags] = xcov(diff(diodeFlipPnts(15:50)), diff(blockFlipPnts(15:50)), 50);
    optimalLag = covarLags(covarValues==max(covarValues));
    if optimalLag < 0;  blockFlipPnts(1:optimalLag*-1)=[];
    elseif optimalLag>0; diodeFlipPnts(1:optimalLag) = [];
    end
    
    initDiff = mean(diodeFlipPnts(25:100) - blockFlipPnts(25:100));
    blockFlipPnts = blockFlipPnts + initDiff;
    blockFlipPnts = blockFlipPnts(knnsearch(blockFlipPnts, diodeFlipPnts));
    blockFlipPnts = blockFlipPnts - initDiff;
end

[covarValues, covarLags] = xcov(diff(diodeFlipPnts), diff(blockFlipPnts), 50);
optimalLag = covarLags(covarValues==max(covarValues));
if optimalLag < 0;  blockFlipPnts(1:optimalLag*-1)=[];
elseif optimalLag>0; diodeFlipPnts(1:optimalLag) = [];
end

minimumFlipNumber = min([length(blockFlipPnts), length(diodeFlipPnts)]);
[blkTConv, stat] = robustfit(blockFlipPnts(2:minimumFlipNumber-2), diodeFlipPnts(2:minimumFlipNumber-2));
if stat.ols_s > 0.02
    warning(['Significant error in fit: ' num2str(stat.ols_s)]);
end
end
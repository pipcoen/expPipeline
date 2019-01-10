function correctedBlock = alignBlock2Timeline(block, timeline, expDef)
%%
switch expDef
    case 'multiSpaceWorld'
        alignType = 'wheel';
        fineTune = {'clicks'; 'flashes'; 'reward'};
end

inputNames = {timeline.hw.inputs.name}';
timelineTime = timeline.rawDAQTimestamps;
sampleRate = 1/diff(timeline.rawDAQTimestamps(1:2));

switch alignType
    case 'wheel'
        smoothWindow = sampleRate/50+1;
        timelinehWeelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
        timelinehWeelPosition(timelinehWeelPosition > 2^31) = timelinehWeelPosition(timelinehWeelPosition > 2^31) - 2^32;
        timelinehWeelPosition = smooth(timelinehWeelPosition,smoothWindow);
        block.inputs.wheelValues = block.inputs.wheelValues-block.inputs.wheelValues(1);
        blockWheelPosition = interp1(block.inputs.wheelTimes, block.inputs.wheelValues, timeline.rawDAQTimestamps, 'linear', 'extrap');
        baseDelay = finddelay(zscore(blockWheelPosition), zscore(timelinehWeelPosition))/sampleRate;
        blockWheelPosition = interp1(block.inputs.wheelTimes+baseDelay, block.inputs.wheelValues, timeline.rawDAQTimestamps, 'linear', 'extrap');
        blockWheelPosition = smooth(blockWheelPosition(:),smoothWindow);
        
        blockWidth = 60*sampleRate;
        sampleCentres = sampleRate*5:sampleRate*5:length(timelineTime);
        blockWheelVelocity = diff(blockWheelPosition);
        timelinehWeelVelocity = diff(timelinehWeelPosition);
        samplePoints = arrayfun(@(x) (x-blockWidth):(x+blockWidth), sampleCentres, 'uni', 0);
        samplePoints = cellfun(@(x) x(x>0 &x<length(timelineTime)), samplePoints, 'uni', 0);
        
        testIdx = cellfun(@(x) sum(abs(blockWheelVelocity(x))), samplePoints)>(5*blockWidth/sampleRate);
        if sum(testIdx) < 50; error('Not enough movment to synchronize using wheel');
        elseif mean(testIdx) < 0.2; warning('Little movement to timeline alignment with wheel will be unreliable');
        end
        samplePoints = samplePoints(testIdx);
        delayValues = cellfun(@(x) finddelay(blockWheelVelocity(x), timelinehWeelVelocity(x), 1000), samplePoints)./sampleRate;
        timelineRefTimes = timelineTime(sampleCentres);
        delayValues = interp1(timelineTime(sampleCentres(testIdx)), delayValues, timelineRefTimes, 'linear', 'extrap');
        delayValues = smooth(delayValues, 0.05, 'rlowess');
        blockRefTimes = movmedian(timelineRefTimes(:)-delayValues-baseDelay, 7)';
        blockRefTimes = interp1(timelineRefTimes(4:end-3), blockRefTimes(4:end-3), timelineRefTimes, 'linear', 'extrap');
    case 'reward'
        blockRefTimes = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
        thresh = max(timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')))/2;
        rewardTrace = timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')) > thresh;
        timelineRefTimes = timeline.rawDAQTimestamps(strfind(rewardTrace', [0 1])+1);
end


fieldList = fieldnames(block.inputs);
fieldList = fieldList(cellfun(@(x) strcmp(x(end-4:end), 'Times'), fieldnames(block.inputs)));
for i = 1:length(fieldList); block.inputs.(fieldList{i}) = interp1(blockRefTimes,timelineRefTimes,block.inputs.(fieldList{i}),'linear','extrap'); end

fieldList = fieldnames(block.outputs);
fieldList = fieldList(cellfun(@(x) strcmp(x(end-4:end), 'Times'), fieldnames(block.outputs)));
for i = 1:length(fieldList); block.outputs.(fieldList{i}) = interp1(blockRefTimes,timelineRefTimes,block.outputs.(fieldList{i}),'linear','extrap'); end

fieldList = fieldnames(block.events);
fieldList = fieldList(cellfun(@(x) strcmp(x(end-4:end), 'Times'), fieldList));
for i = 1:length(fieldList); block.events.(fieldList{i}) = interp1(blockRefTimes,timelineRefTimes,block.events.(fieldList{i}),'linear','extrap'); end

cutOff = round(block.events.newTrialTimes(1)*sampleRate-(0.05*sampleRate));
timeline.rawDAQData(1:cutOff,:) = repmat(timeline.rawDAQData(cutOff,:), cutOff, 1);

%%
if contains('reward', fineTune)
    blockRewardTimes = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
    thresh = max(timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')))/2;
    rewardTrace = timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')) > thresh;
    timelineRewardTimes = timeline.rawDAQTimestamps(strfind(rewardTrace', [0 1])+1);
    if length(timelineRewardTimes)~=length(blockRewardTimes); error('There should always be an equal number of reward signals'); end
    block.outputs.rewardTimes(block.outputs.rewardValues > 0) = timelineRewardTimes(prc.nearestPoint(blockRewardTimes,timelineRewardTimes));
end

if contains('clicks', fineTune)
    if ~contains('audioOut', inputNames); warning('No audio output in timeline so cannot fine tune'); end
    timelineClickTrace = [0;diff(detrend(timeline.rawDAQData(:,strcmp(inputNames, 'audioOut'))))];
    [~, thresh] = kmeans(timelineClickTrace,3);
    timelineClickOn = timelineTime(strfind((timelineClickTrace>max(thresh)*0.5)', [0 1]));
    timelineClickOff = timelineTime(strfind((timelineClickTrace<min(thresh)*0.5)', [0 1]));
    detectedDuration = round(mean(timelineClickOff-timelineClickOn)*1000);
    
    if length(timelineClickOn)~=length(timelineClickOff); error('There should always be an equal number on/off signals for clicks'); end
    if detectedDuration~=unique([block.paramsValues.clickDuration])*1000; error('Diff in detected and requested click durations'); end
    
    audStimOnOffTimesValues = sortrows([[timelineClickOn';timelineClickOff'] [timelineClickOn'*0+1; timelineClickOff'*0]],1);    
    block.events.audStimOnOffTimes = audStimOnOffTimesValues(:,1)';
    block.events.audStimOnOffValues = audStimOnOffTimesValues(:,2)';
    
    largeAudGaps = sort([find(diff([0; audStimOnOffTimesValues(:,1)])>0.5); find(diff([audStimOnOffTimesValues(:,1); 10e10])>0.5)]);
    block.events.audStimPeriodOnOffTimes = audStimOnOffTimesValues(largeAudGaps,1)';
    block.events.audStimPeriodOnOffValues = audStimOnOffTimesValues(largeAudGaps,2)';
    
    [compareIndex] = prc.nearestPoint(block.events.stimPeriodOnOffTimes, block.events.audStimPeriodOnOffTimes);
    if any(compareIndex-(1:numel(compareIndex))); fprintf('Error in matching visual stimulus start and end times \n'); keyboard; end
end

%%
if contains('flashes', fineTune)
    % Change visual stimulus times to timeline version
    photoDiodeTrace = timeline.rawDAQData(:,strcmp(inputNames, 'photoDiode'));
    [~, thresh] = kmeans(photoDiodeTrace,2);
    thresh = [min(thresh) + range(thresh)*0.25;  max(thresh) - range(thresh)*0.25];
    photoDiodeFlipOn = sort([strfind(photoDiodeTrace'>thresh(1), [0 1]), strfind(photoDiodeTrace'>thresh(2), [0 1])]);
    photoDiodeFlipOff = sort([strfind(photoDiodeTrace'<thresh(1), [0 1]), strfind(photoDiodeTrace'<thresh(2), [0 1])]);
    photoDiodeFlips = sort([photoDiodeFlipOn photoDiodeFlipOff]);
    photoDiodeFlips([strfind(ismember(photoDiodeFlips, photoDiodeFlipOn), [1 1])+1 strfind(ismember(photoDiodeFlips, photoDiodeFlipOff), [1 1])+1]) = [];
    photoDiodeFlipTimes = timeline.rawDAQTimestamps(photoDiodeFlips)';
    %%
    largeVisGaps = photoDiodeFlipTimes(sort([find(diff([0; photoDiodeFlipTimes])>0.5); find(diff([photoDiodeFlipTimes; 10e10])>0.5)]));
    block.events.visStimPeriodOnOffValues = block.events.stimPeriodOnOffValues;
    block.events.visStimPeriodOnOffTimes = largeVisGaps(1:length(block.events.stimPeriodOnOffValues))';
    
    [compareIndex] = prc.nearestPoint(block.events.stimPeriodOnOffTimes, block.events.visStimPeriodOnOffTimes);
    if any(compareIndex-(1:numel(compareIndex))); fprintf('Error in matching visual stimulus start and end times \n'); keyboard; end
end
correctedBlock = block;
correctedBlock.fineTuned = 1;
end
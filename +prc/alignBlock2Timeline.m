function [correctedBlock, aligned] = alignBlock2Timeline(block, timeline, expDef)
%%
inputNames = {timeline.hw.inputs.name}';
timelineTime = timeline.rawDAQTimestamps;
sampleRate = 1/diff(timeline.rawDAQTimestamps(1:2));


switch expDef
    case 'multiSpaceWorld'
        alignType = 'wheel';
        fineTune = {'clicksfine'; 'flashes'; 'reward'};
    case 'multiSpaceWorldPassive'
        alignType = 'wheel';
        fineTune = {'clicksfine'; 'flashesfine'; 'reward'};
end


switch alignType
    case 'wheel'
        smoothWindow = sampleRate/10+1;
        timelinehWeelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
        timelinehWeelPosition(timelinehWeelPosition > 2^31) = timelinehWeelPosition(timelinehWeelPosition > 2^31) - 2^32;
        timelinehWeelPosition = smooth(timelinehWeelPosition,smoothWindow);
        block.inputs.wheelValues = block.inputs.wheelValues-block.inputs.wheelValues(1);
        blockWheelPosition = interp1(block.inputs.wheelTimes, block.inputs.wheelValues, timeline.rawDAQTimestamps, 'linear', 'extrap');
        blockWheelPosition = smooth(blockWheelPosition,smoothWindow);
        baseDelay = finddelay(diff(blockWheelPosition), diff(timelinehWeelPosition))/sampleRate;
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

trialStEnTimes = [block.events.newTrialTimes(1:length(block.events.endTrialTimes))', block.events.endTrialTimes'];
stimStartBlock = block.events.stimPeriodOnOffTimes(block.events.stimPeriodOnOffValues==1);
%%
 if contains('reward', fineTune)
    blockRewardTimes = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
    thresh = max(timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')))/2;
    rewardTrace = timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')) > thresh;
    timelineRewardTimes = timeline.rawDAQTimestamps(strfind(rewardTrace', [0 1])+1);
    if length(timelineRewardTimes)~=length(blockRewardTimes); error('There should always be an equal number of reward signals'); end
    block.outputs.rewardTimes(block.outputs.rewardValues > 0) = timelineRewardTimes(prc.nearestPoint(blockRewardTimes,timelineRewardTimes));
    aligned.rewardTimes = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
    aligned.rewardTimes = aligned.rewardTimes(:);
end

if any(contains(fineTune, 'clicks')) && contains('audioOut', inputNames)
    timelineClickTrace = [0;diff(detrend(timeline.rawDAQData(:,strcmp(inputNames, 'audioOut'))))];
    [~, thresh] = kmeans(timelineClickTrace,3);
    timelineClickOn = timelineTime(strfind((timelineClickTrace>max(thresh)*0.5)', [0 1]));
    timelineClickOff = timelineTime(strfind((timelineClickTrace<min(thresh)*0.5)', [0 1]));
    detectedDuration = round(mean(timelineClickOff-timelineClickOn)*1000);
    
    if length(timelineClickOn)~=length(timelineClickOff); error('There should always be an equal number on/off signals for clicks'); end
    if abs(detectedDuration-(unique([block.paramsValues.clickDuration])*1000))>3; error('Diff in detected and requested click durations'); end
    
    aStimOnOffTV = sortrows([[timelineClickOn';timelineClickOff'] [timelineClickOn'*0+1; timelineClickOff'*0]],1);    
    block.events.audStimOnOffTimes = aStimOnOffTV(:,1)';
    block.events.audStimOnOffValues = aStimOnOffTV(:,2)';
    aligned.audStimOnOff = [aStimOnOffTV(aStimOnOffTV(:,2)==1,1) aStimOnOffTV(aStimOnOffTV(:,2)==0,1)];
    
    largeAudGaps = sort([find(diff([0; aStimOnOffTV(:,1)])>0.5); find(diff([aStimOnOffTV(:,1); 10e10])>0.5)]);
    block.events.audStimPeriodOnOffTimes = aStimOnOffTV(largeAudGaps,1)';
    block.events.audStimPeriodOnOffValues = aStimOnOffTV(largeAudGaps,2)';
    
    %% Sanity check (should be match between stim starts from block and from timeline)
    audstimStartTimeline = block.events.audStimPeriodOnOffTimes(block.events.audStimPeriodOnOffValues==1);
    nonAudTrials = [block.paramsValues.audAmplitude]==0;
    [compareIndex] = prc.nearestPoint(stimStartBlock(~nonAudTrials), audstimStartTimeline);
    if any(compareIndex-(1:numel(compareIndex))); fprintf('Error in matching visual stimulus start and end times \n'); keyboard; end
    
    aStimOnOffTV = aStimOnOffTV(largeAudGaps,:);
    aligned.audStimPeriodOnOff = [aStimOnOffTV(aStimOnOffTV(:,2)==1,1) aStimOnOffTV(aStimOnOffTV(:,2)==0,1)];
elseif contains('clicks', fineTune)
    warning('No audio output in timeline so cannot fine tune');
end

%%
if any(contains(fineTune, 'flashes'))
    % Change visual stimulus times to timeline version
    photoDiodeTrace = timeline.rawDAQData(:,strcmp(inputNames, 'photoDiode'));
    [~, thresh] = kmeans(photoDiodeTrace,2);
    thresh = [min(thresh) + range(thresh)*0.25;  max(thresh) - range(thresh)*0.25];
    photoDiodeFlipOn = sort([strfind(photoDiodeTrace'>thresh(1), [0 1]), strfind(photoDiodeTrace'>thresh(2), [0 1])]);
    photoDiodeFlipOff = sort([strfind(photoDiodeTrace'<thresh(1), [0 1]), strfind(photoDiodeTrace'<thresh(2), [0 1])]);
    photoDiodeFlips = sort([photoDiodeFlipOn photoDiodeFlipOff]);
    photoDiodeFlips([strfind(ismember(photoDiodeFlips, photoDiodeFlipOn), [1 1])+1 strfind(ismember(photoDiodeFlips, photoDiodeFlipOff), [1 1])+1]) = [];
    photoDiodeFlipTimes = timeline.rawDAQTimestamps(photoDiodeFlips)';
    photoDiodeFlipTimes(find(diff(photoDiodeFlipTimes)<(12/1000))+1) = [];
    
    largeVisGaps = photoDiodeFlipTimes(sort([find(diff([0; photoDiodeFlipTimes])>0.5); find(diff([photoDiodeFlipTimes; 10e10])>0.5)]));
    zeroContrastTrials = arrayfun(@(x) max(abs(x.visContrast)),block.paramsValues)'==0;
    largeGapsByTrial = arrayfun(@(x,y) find(largeVisGaps>x & largeVisGaps<y), trialStEnTimes(:,1), trialStEnTimes(:,2), 'uni', 0);
    largeVisGaps(cell2mat(largeGapsByTrial(zeroContrastTrials))) = [];
    
    %% Sanity check (should be match between stim starts from block and from timeline)
    [compareIndex] = prc.nearestPoint(stimStartBlock(~zeroContrastTrials), largeVisGaps(1:2:end)');
    if any(compareIndex-(1:numel(compareIndex)))
        fprintf('WARNING: problem matching visual stimulus start and end times \n');
        fprintf('Will try removing points that do not match stimulus starts \n');
        
        [~, nearestPoint] = prc.nearestPoint(largeVisGaps, block.events.stimPeriodOnOffTimes);
        largeVisGaps(nearestPoint>0.5) = [];
        
        [compareIndex] = prc.nearestPoint(block.events.stimPeriodOnOffTimes, largeVisGaps(1:length(block.events.stimPeriodOnOffValues))');
        if any(compareIndex-(1:numel(compareIndex))); fprintf('Error in matching visual stimulus start and end times \n'); keyboard; end
        block.events.visStimPeriodOnOffTimes = largeVisGaps(1:length(block.events.stimPeriodOnOffValues))';
    end
    %%
    block.events.visStimPeriodOnOffValues = 0*largeVisGaps'+1;
    block.events.visStimPeriodOnOffValues(2:2:end) = 0;
    block.events.visStimPeriodOnOffTimes = largeVisGaps';

    vStimOnOffTV = [block.events.visStimPeriodOnOffTimes' block.events.visStimPeriodOnOffValues'];
    vStimOnOffTV = vStimOnOffTV(1:find(vStimOnOffTV(:,2)==0,1,'last'),:);
    aligned.visStimPeriodOnOff = [vStimOnOffTV(vStimOnOffTV(:,2)==1,1) vStimOnOffTV(vStimOnOffTV(:,2)==0,1)];
end

if any(contains(fineTune, 'flashesfine'))
    photoFlipsByTrial = arrayfun(@(x,y) find(photoDiodeFlipTimes>=x & photoDiodeFlipTimes<=y), aligned.visStimPeriodOnOff(:,1), aligned.visStimPeriodOnOff(:,2), 'uni', 0);
    expectedFlashTrainLength = [block.paramsValues.clickRate]'.*[block.paramsValues.responseWindow]'*2;
    if any(expectedFlashTrainLength(~zeroContrastTrials)-cellfun(@length,photoFlipsByTrial)); fprintf('Error in flashesfine tuning \n'); keyboard; end
    
    block.events.visStimOnOffTimes = sort(photoDiodeFlipTimes(cell2mat(photoFlipsByTrial)))';
    block.events.visStimOnOffValues = photoDiodeFlipTimes(cell2mat(photoFlipsByTrial))'*0+1;
    block.events.visStimOnOffValues(2:2:end) = 0;
    vStimOnOffTV = [block.events.visStimOnOffTimes' block.events.visStimOnOffValues'];
    aligned.visStimOnOff = [vStimOnOffTV(vStimOnOffTV(:,2)==1,1) vStimOnOffTV(vStimOnOffTV(:,2)==0,1)];
end

correctedBlock = block;
correctedBlock.fineTuned = 1;
end
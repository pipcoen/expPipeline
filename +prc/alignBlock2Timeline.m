function [correctedBlock, aligned] = alignBlock2Timeline(block, timeline, expDef)
%% A funciton align the block times to the timeline times for ephys, fusi, etc. experiments

%INPUTS(default values)
%block(required)----------The block file to be aligned
%timeline(required)-------The timeline file to be aligned
%expDef(required)---------A string that defines the expDef used in the experiment

%OUTPUTS
%correctedBlock-----------The corrected block, with times now matching timeline
%aligned------------------
%"eph" is a structure with ephys-specific information
%"whoD"is a cell will the names of all variables in the file (loading this is quicker than asking matlab to check)

%% Get a set of corresponding timeline and block timepoints for interpolation
inputNames = {timeline.hw.inputs.name}';                      %List of inputs to the timeline file
timelineTime = timeline.rawDAQTimestamps;                     %Timestamps in the timeline file
sR = 1/diff(timeline.rawDAQTimestamps(1:2));          %Timeline sample rate

%Establish whether there is a record in timeline of the audio output
if contains('audioOut', inputNames); audInput = 'audioOut';
else, audInput = 0;
end

%Which features are extracted from timline depends on the expDef. Individual on/off times for each click are detected only in Passive.
%"alignType" defines the type of signal to use in the initial alignment. It can be "wheel", "reward" or "photoDiode." "wheel" is best for active
%experiments (because it occurs most often) and photoDiode is best for passive experiments (because wheel movement is likely minimal).
%"fineTune" defines the timeline features to extract from each expDef.
switch expDef
    case 'multiSpaceWorld'
        alignType = 'wheel';
        fineTune = {'clicksfine'; 'flashes'; 'reward'; 'movements'; 'wheelTraceTimeValue'};
        trialGapThresh = 1;
    case {'multiSpaceWorldPassive'; 'multiSpaceWorldPassiveTwoSpeakers'}
        if contains('photoDiode', inputNames); alignType = 'photoDiode';
        elseif contains('rotaryEncoder', inputNames); alignType = 'wheel';
        else, error('What am I supposed to use to align block with timeline?!')
        end
        trialGapThresh = 1./max([block.paramsValues.clickRate])*3;
        fineTune = {'clicksfine'; 'flashesfine'; 'reward'; 'wheelTraceTimeValue'};
end

switch alignType
    case 'wheel' %Get interpolation points using the wheel data
        %Unrap the wheel trace (it is circular) and then smooth it. Smoothing is important because covariance will not work otherwise
        smthWin = sR/10+1;
        timelinehWeelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
        timelinehWeelPosition(timelinehWeelPosition > 2^31) = timelinehWeelPosition(timelinehWeelPosition > 2^31) - 2^32;
        timelinehWeelPositionSmooth = smooth(timelinehWeelPosition,smthWin);

        %Make initial block time zero (as with timeline) and then interpolate with timeline times and smooth to give same number of points etc.
        block.inputs.wheelValues = block.inputs.wheelValues-block.inputs.wheelValues(1);
        blockWheelPosition = interp1(block.inputs.wheelTimes, block.inputs.wheelValues, timeline.rawDAQTimestamps, 'linear', 'extrap');
        blockWheelPosition = smooth(blockWheelPosition,smthWin);
        
        %Find the overall delay with the entire trace. Then reinterpolate, accounting for this delay.
        baseDelay = finddelay(diff(blockWheelPosition), diff(timelinehWeelPositionSmooth))/sR;
        blockWheelPosition = interp1(block.inputs.wheelTimes+baseDelay, block.inputs.wheelValues, timeline.rawDAQTimestamps, 'linear', 'extrap');
        blockWheelPosition = smooth(blockWheelPosition(:),smthWin);
        
        %Get the vectors for 20s sections of the blockWheelVelocity and timelinehWeelVelocity
        blockWidth = 20*sR;
        sampleCentres = sR*5:sR*10:length(timelineTime);
        blockWheelVelocity = diff(blockWheelPosition);
        timelinehWeelVelocity = diff(timelinehWeelPosition);
        samplePoints = arrayfun(@(x) (x-blockWidth):(x+blockWidth), sampleCentres, 'uni', 0);
        samplePoints = cellfun(@(x) x(x>0 & x<length(timelineTime)), samplePoints, 'uni', 0);
        
        %Check that there is enough wheel movement to make the alignement (based on absolute velocity)
        testIdx = cellfun(@(x) sum(abs(blockWheelVelocity(x))), samplePoints)>(5*blockWidth/sR);
        if mean(testIdx) < 0.2; error('Not enough movment to synchronize using wheel');
        elseif mean(testIdx) < 0.2; warning('Little movement so timeline alignment with wheel will be unreliable');
        end
        
        %Go through each subsection and detect the offset between block and timline
        samplePoints = samplePoints(testIdx);
        delayValues = cellfun(@(x) finddelay(blockWheelVelocity(x), timelinehWeelVelocity(x), 1000), samplePoints)./sR;
        
        %Use a smoothed median to select the evolving delat values, and use these to calculate the evolving reference points for block and timeline
        timelineRefTimes = timelineTime(sampleCentres);
        delayValues = interp1(timelineTime(sampleCentres(testIdx)), delayValues, timelineRefTimes, 'linear', 'extrap');
        delayValues = smooth(delayValues, 0.05, 'rlowess');
        blockRefTimes = movmedian(timelineRefTimes(:)-delayValues-baseDelay, 7)';
        blockRefTimes = interp1(timelineRefTimes(4:end-3), blockRefTimes(4:end-3), timelineRefTimes, 'linear', 'extrap');
        block.alignment = 'wheel';
        
    case 'reward' %Get interpolation points using the reward data
        %Rewards are very obvious, but infrequent and sometimes completely absent. But if there is no other option, simply detect the rewards in
        %timesline, and use the rewardTimes from the block outputs to get reference points
        blockRefTimes = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
        thresh = max(timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')))/2;
        rewardTrace = timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho')) > thresh;
        timelineRefTimes = timeline.rawDAQTimestamps(strfind(rewardTrace', [0 1])+1);
        if length(timelineRefTimes)>length(blockRefTimes); timelineRefTimes = timelineRefTimes(2:end); end
        block.alignment = 'reward';
        
    case 'photoDiode' %Get interpolation points using the photodiode data
        %Note, the photodiode *should* vary between two values, but it often doesn't. For reasons we don't understand, it sometimes goes to grey, and
        %sometimes you skip changes that are recorded in the block etc. This is another reason I prefer to use the wheel. But this method has proved
        %reasonably reliable if the wheel isn't suitable.
        
        %Extract photodiode trace and get repeated values by using kmeans. Get the lower and upper thersholds from this range.
        photoDiodeTrace = timeline.rawDAQData(:,strcmp(inputNames, 'photoDiode'));
        [~, thresh] = kmeans(photoDiodeTrace,5);
        thresh = [min(thresh) + range(thresh)*0.2;  max(thresh) - range(thresh)*0.2];
        
        %Find flips based on these thresholds.
        photoDiodeFlipOn = sort([strfind(photoDiodeTrace'>thresh(1), [0 1]), strfind(photoDiodeTrace'>thresh(2), [0 1])]);
        photoDiodeFlipOff = sort([strfind(photoDiodeTrace'<thresh(1), [0 1]), strfind(photoDiodeTrace'<thresh(2), [0 1])]);
        photoDiodeFlips = sort([photoDiodeFlipOn photoDiodeFlipOff]);
        
        %Remove cases where two flips in the same direction appear in succession (you can't flip to white twice in a row)
        photoDiodeFlips([strfind(ismember(photoDiodeFlips, photoDiodeFlipOn), [1 1])+1 strfind(ismember(photoDiodeFlips, photoDiodeFlipOff), [1 1])+1]) = [];
        
        %Get corresponding flip times. Remove any that would be faster than 60Hz (screen refresh rate)
        photoDiodeFlipTimes = timeline.rawDAQTimestamps(photoDiodeFlips)';
        photoDiodeFlipTimes(find(diff(photoDiodeFlipTimes)<(12/1000))+1) = [];
        blockRefTimes = block.stimWindowUpdateTimes(diff(block.stimWindowUpdateTimes)>0.49);
        timelineRefTimes = photoDiodeFlipTimes(diff(photoDiodeFlipTimes)>0.49);
        
        %Use "prc.try2alignVectors" to deal with cases where the timeline and block flip times are different lengths, or have large differences. I
        %have found this to solve all problems like this. However, I have also found it to be critical (the photodiode is just too messy otherwise)
        if length(blockRefTimes) ~= length(timelineRefTimes)
            [timelineRefTimes, blockRefTimes] = prc.try2alignVectors(timelineRefTimes, blockRefTimes, 0.25);
        elseif any(abs((blockRefTimes-blockRefTimes(1)) - (timelineRefTimes-timelineRefTimes(1)))>0.5)
            [timelineRefTimes, blockRefTimes] = prc.try2alignVectors(timelineRefTimes, blockRefTimes, 0.25);
        end
        block.alignment = 'photodiode';
        if length(blockRefTimes) ~= length(timelineRefTimes)
            error('Photodiode alignment error');
        end
end

%% Use referece points to interpolate various fields of the block file
%This does exactly the same thing for "inputs", "outputs" and "events" based on the field ending in "Times"
fieldList = fieldnames(block.inputs);
fieldList = fieldList(cellfun(@(x) strcmp(x(end-4:end), 'Times'), fieldnames(block.inputs)));
for i = 1:length(fieldList); block.inputs.(fieldList{i}) = interp1(blockRefTimes,timelineRefTimes,block.inputs.(fieldList{i}),'linear','extrap'); end

fieldList = fieldnames(block.outputs);
fieldList = fieldList(cellfun(@(x) strcmp(x(end-4:end), 'Times'), fieldnames(block.outputs)));
for i = 1:length(fieldList); block.outputs.(fieldList{i}) = interp1(blockRefTimes,timelineRefTimes,block.outputs.(fieldList{i}),'linear','extrap'); end

fieldList = fieldnames(block.events);
fieldList = fieldList(cellfun(@(x) strcmp(x(end-4:end), 'Times'), fieldList));
for i = 1:length(fieldList); block.events.(fieldList{i}) = interp1(blockRefTimes,timelineRefTimes,block.events.(fieldList{i}),'linear','extrap'); end

%Flatten timeline traces prior to the first trial starting (avoids errneous event idenitification
cutOff = round(block.events.newTrialTimes(1)*sR-(0.05*sR));
timeline.rawDAQData(1:cutOff,:) = repmat(timeline.rawDAQData(cutOff,:), cutOff, 1);

%% Prepare to extract events from timeline
%Get the start/end times of the trials, and the stimulus start times from the blocks.
trialStEnTimes = [block.events.newTrialTimes(1:length(block.events.endTrialTimes))', block.events.endTrialTimes'];
stimStartBlock = block.events.stimPeriodOnOffTimes(block.events.stimPeriodOnOffValues==1);

%% Extract the reward times from timeline. This is straightforward as the reward echo is very reliable. Just thershold and detect peaks.
if contains('reward', fineTune) && contains('rewardEcho',inputNames)
    blockRewardTimes = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
    rewardTrace = mat2gray(timeline.rawDAQData(:,strcmp(inputNames, 'rewardEcho'))) > 0.5;
    timelineRewardTimes = timeline.rawDAQTimestamps(strfind(rewardTrace', [0 1])+1);
    if length(timelineRewardTimes)~=length(blockRewardTimes); error('There should always be an equal number of reward signals'); end
    block.outputs.rewardTimes(block.outputs.rewardValues > 0) = timelineRewardTimes(prc.nearestPoint(blockRewardTimes,timelineRewardTimes));
    aligned.rewardTimes = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
    aligned.rewardTimes = aligned.rewardTimes(:);
end
if contains('reward', fineTune) && ~contains('rewardEcho',inputNames); warning('No reward input echo... skipping'); end

%% Extract audio clicks (these are pretty reliable, so can extract every click)
if any(contains(fineTune, 'clicks')) && ischar(audInput)
    %Detrend timeline trace, threshold using kmeans, detect onsets and offsets of sound, estimate duration from this.
    timelineClickTrace = [0;diff(detrend(timeline.rawDAQData(:,strcmp(inputNames, audInput))))];
    [~, thresh] = kmeans(timelineClickTrace,5);
    timelineClickOn = timelineTime(strfind((timelineClickTrace>max(thresh)*0.25)', [0 1]));
    timelineClickOff = timelineTime(strfind((timelineClickTrace<min(thresh)*0.25)', [0 1]));
    detectedDuration = round(mean(timelineClickOff-timelineClickOn)*1000);
    
    %Sanity check: same number of onsets and offsets, check that detected duration matches the duration parameter (assumed constant here)
    if length(timelineClickOn)~=length(timelineClickOff); error('There should always be an equal number on/off signals for clicks'); end
    if abs(detectedDuration-(unique([block.paramsValues.clickDuration])*1000))>3; error('Diff in detected and requested click durations'); end
    
    %Create vector that is sorted by time: [onset time, offset time, 1, 0] and find large gaps between successive onsets (stimulus period onsets)
    aStimOnOffTV = sortrows([[timelineClickOn';timelineClickOff'] [timelineClickOn'*0+1; timelineClickOff'*0]],1);  
    largeAudGaps = sort([find(diff([0; aStimOnOffTV(:,1)])>trialGapThresh); find(diff([aStimOnOffTV(:,1); 10e10])>trialGapThresh)]);

    %%%%%%%%%%%%%%%%%%%%STILL NEEDS TO BE COMMENTED BELOW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sanity check (should be match between stim starts from block and from timeline)
    audstimStartTimeline = aStimOnOffTV(largeAudGaps,1);
    audstimStartTimeline = audstimStartTimeline(aStimOnOffTV(largeAudGaps,2)==1);
    nonAudTrials = [block.paramsValues.audAmplitude]==0; nonAudTrials = nonAudTrials(1:length(stimStartBlock));
    [compareIndex] = prc.nearestPoint(stimStartBlock(~nonAudTrials), audstimStartTimeline);
    audError = 0;
    if any(compareIndex-(1:numel(compareIndex)))
        audError = 1;
        [compareIndex] = prc.nearestPoint(stimStartBlock(~nonAudTrials), audstimStartTimeline(~nonAudTrials));
        if ~any(compareIndex-(1:numel(compareIndex))) && length(largeAudGaps)/2 == length(nonAudTrials)
            fprintf('WARNING: Detected that AmpAud = 0 trials still generated signals in Timeline. Will remove these \n')
            
            largeAudGaps([find(nonAudTrials)*2-1 find(nonAudTrials)*2]) = [];
            audstimStartTimeline = aStimOnOffTV(largeAudGaps,1);
            audstimStartTimeline = audstimStartTimeline(aStimOnOffTV(largeAudGaps,2)==1);
            [compareIndex] = prc.nearestPoint(stimStartBlock(~nonAudTrials), audstimStartTimeline);
            if ~any(compareIndex-(1:numel(compareIndex)))
                audError = 0; 
                keepIdx = cell2mat(arrayfun(@(x) largeAudGaps(x):largeAudGaps(x+1), 1:2:length(largeAudGaps), 'uni', 0));
                aStimOnOffTV = aStimOnOffTV(keepIdx,:);
                largeAudGaps = sort([find(diff([0; aStimOnOffTV(:,1)])>trialGapThresh); find(diff([aStimOnOffTV(:,1); 10e10])>trialGapThresh)]);
            end
        end
    end
    if audError; fprintf('Error in matching auditory stimulus start and end times \n'); keyboard; end
    
    block.events.audStimOnOffTimes = aStimOnOffTV(:,1)';
    block.events.audStimOnOffValues = aStimOnOffTV(:,2)';
    aligned.audStimOnOff = [aStimOnOffTV(aStimOnOffTV(:,2)==1,1) aStimOnOffTV(aStimOnOffTV(:,2)==0,1)];
    block.events.audStimPeriodOnOffTimes = aStimOnOffTV(largeAudGaps,1)';
    block.events.audStimPeriodOnOffValues = aStimOnOffTV(largeAudGaps,2)';
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
    
    trialGapThresh = 1./max([block.paramsValues.clickRate])*3;
    largeVisGaps = photoDiodeFlipTimes(sort([find(diff([0; photoDiodeFlipTimes])>trialGapThresh); find(diff([photoDiodeFlipTimes; 10e10])>trialGapThresh)]));
    zeroContrastTrials = arrayfun(@(x) max(abs(x.visContrast)),block.paramsValues)'==0;
    zeroContrastTrials = zeroContrastTrials(1:size(trialStEnTimes,1));
    largeGapsByTrial = arrayfun(@(x,y) largeVisGaps(largeVisGaps>x & largeVisGaps<y), trialStEnTimes(:,1), trialStEnTimes(:,2), 'uni', 0);
    largeVisGaps = cell2mat(largeGapsByTrial(~zeroContrastTrials));
    largeVisGaps = [largeVisGaps(1:2:end) largeVisGaps(2:2:end)];
    
    %% Sanity check (should be match between stim starts from block and from timeline)
    [compareIndex] = prc.nearestPoint(stimStartBlock(~zeroContrastTrials), largeVisGaps(:,1)');
    if any(compareIndex-(1:numel(compareIndex)))
        fprintf('WARNING: problem matching visual stimulus start and end times \n');
        fprintf('Will try removing points that do not match stimulus starts \n');
        
        [~, nearestPoint] = prc.nearestPoint(largeVisGaps(:,1), stimStartBlock(~zeroContrastTrials));
        largeVisGaps(nearestPoint>0.5,:) = [];
        
        [compareIndex] = prc.nearestPoint(stimStartBlock(~zeroContrastTrials), largeVisGaps(:,1)')';
        if any(compareIndex-(1:numel(compareIndex))'); fprintf('Error in matching visual stimulus start and end times \n'); keyboard; end
    end
    %%
    block.events.visStimPeriodOnOffValues = 0*sort(largeVisGaps(:))'+1;
    block.events.visStimPeriodOnOffValues(2:2:end) = 0;
    block.events.visStimPeriodOnOffTimes = sort(largeVisGaps(:))';
    
    vStimOnOffTV = [block.events.visStimPeriodOnOffTimes' block.events.visStimPeriodOnOffValues'];
    vStimOnOffTV = vStimOnOffTV(1:find(vStimOnOffTV(:,2)==0,1,'last'),:);
    aligned.visStimPeriodOnOff = [vStimOnOffTV(vStimOnOffTV(:,2)==1,1) vStimOnOffTV(vStimOnOffTV(:,2)==0,1)];
end

if any(contains(fineTune, 'flashesfine'))
    photoFlipsByTrial = arrayfun(@(x,y) find(photoDiodeFlipTimes>=x & photoDiodeFlipTimes<=y), aligned.visStimPeriodOnOff(:,1), aligned.visStimPeriodOnOff(:,2), 'uni', 0);
    expectedFlashTrainLength = [block.paramsValues.clickRate]'.*[block.paramsValues.responseWindow]'*2;
    misMatchFlashtrain = expectedFlashTrainLength(~zeroContrastTrials)-cellfun(@length,photoFlipsByTrial);
    if any(misMatchFlashtrain)
        photoFlipsByTrial(misMatchFlashtrain~=0) = [];
        fprintf('Warning: Removing flash times for trials that do not match predicted flash length \n');
    end
    
    block.events.visStimOnOffTimes = sort(photoDiodeFlipTimes(cell2mat(photoFlipsByTrial)))';
    block.events.visStimOnOffValues = photoDiodeFlipTimes(cell2mat(photoFlipsByTrial))'*0+1;
    block.events.visStimOnOffValues(2:2:end) = 0;
    vStimOnOffTV = [block.events.visStimOnOffTimes' block.events.visStimOnOffValues'];
    aligned.visStimOnOff = [vStimOnOffTV(vStimOnOffTV(:,2)==1,1) vStimOnOffTV(vStimOnOffTV(:,2)==0,1)];
end

if any(contains(fineTune, 'movements'))
    responseMadeIdx = block.events.feedbackValues(1:size(trialStEnTimes,1)) ~= 0;
    
    timelineVisOnset = prc.indexByTrial(trialStEnTimes, aligned.visStimPeriodOnOff(:,1), aligned.visStimPeriodOnOff(:,1));
    timelineVisOnset(cellfun(@isempty, timelineVisOnset)) = deal({nan});
    timelineAudOnset = prc.indexByTrial(trialStEnTimes, aligned.audStimPeriodOnOff(:,1), aligned.audStimPeriodOnOff(:,1));
    timelineAudOnset(cellfun(@isempty, timelineAudOnset)) = deal({nan});
    timelineStimOnset = nanmin(cell2mat([timelineVisOnset timelineAudOnset]), [],2);
    timelineStimOnsetIdx = round(timelineStimOnset(responseMadeIdx)*sR)-100;  %We want to know if movements were initiatied 100ms preceding stim onset
    %%
    stimOnsetIdx = timelineStimOnsetIdx;
    closedLoopPeriodIdx = round(block.events.closedLoopOnOffTimes*sR);
    closedLoopValues = block.events.closedLoopOnOffValues(1:find(block.events.closedLoopOnOffValues==0, 1, 'last'));
    
    if ~exist('timelinehWeelPosition', 'var')
        timelinehWeelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
        timelinehWeelPosition(timelinehWeelPosition > 2^31) = timelinehWeelPosition(timelinehWeelPosition > 2^31) - 2^32;
    end
    wheel = timelinehWeelPosition;
    move4Response = wheel(closedLoopPeriodIdx(closedLoopValues==0)) - wheel(closedLoopPeriodIdx(closedLoopValues==1));
    wheelThresh = median(abs(move4Response(responseMadeIdx)))*0.4;
    threshCrsIdx = arrayfun(@(x,y) max([nan find(abs(wheel(x:(x+(sR*2)))-wheel(x))>wheelThresh,1)+x]), stimOnsetIdx);

    sumWin = 51;
    velThresh  = sR*3/sumWin; %Note: This is somewhat arbitrary. Would be different for different rigs...
    wheelVel = diff([0; wheel])*sR;
    posVelScan = conv(wheelVel.*double(wheelVel>0) - double(wheelVel<0)*1e6, [ones(1,sumWin) zeros(1,sumWin-1)]./sumWin, 'same').*(wheelVel~=0);
    negVelScan = conv(wheelVel.*double(wheelVel<0) + double(wheelVel>0)*1e6, [ones(1,sumWin) zeros(1,sumWin-1)]./sumWin, 'same').*(wheelVel~=0);
    
    moveOnsetIdx = [strfind((posVelScan'>=velThresh), [0,1]) -1*strfind((-1*negVelScan'>=velThresh), [0,1])];
    movingIdx = sort([find(posVelScan'>=velThresh) find((-1*negVelScan'>=velThresh))]);
    [~, srtIdx] = sort(abs(moveOnsetIdx));
    moveOnsetSign = sign(moveOnsetIdx(srtIdx));
    moveOnsetIdx = abs(moveOnsetIdx(srtIdx));
    
    
    %"firstMoveTimes" are the first onsets occuring after stimOnsetIdx. "largeMoveTimes" are the first onsets occuring after stimOnsetIdx that match the
    %sign of the threshold crossing defined earlier. Eliminate any that are longer than 1.5s, as these would be timeouts. Also, remove onsets when the
    %mouse was aready moving at the time of the stimulus onset (impossible to get an accurate movement onset time in this case)
    firstMoveIdx =  arrayfun(@(x) max([nan moveOnsetIdx(find(moveOnsetIdx>x, 1))]), stimOnsetIdx);
    firstMoveIdx((firstMoveIdx-stimOnsetIdx)>1.5*sR) = nan;
    firstMoveIdx(ismember(stimOnsetIdx, movingIdx)) = nan;
    firstMoveDir = arrayfun(@(x) max([nan moveOnsetSign(find(moveOnsetIdx>x, 1))]), stimOnsetIdx);
    firstMoveDir(isnan(firstMoveIdx)) = nan;
    firstMoveDir = ((firstMoveDir==-1)+1).*(abs(firstMoveDir));
    
    preThreshOnsets = arrayfun(@(x,y) moveOnsetIdx(moveOnsetIdx>=x & moveOnsetIdx<=y), stimOnsetIdx, threshCrsIdx, 'uni', 0);
    preThreshDirs = arrayfun(@(x,y) moveOnsetSign(moveOnsetIdx>=x & moveOnsetIdx<=y), stimOnsetIdx, threshCrsIdx, 'uni', 0);
    reliableMoves = firstMoveDir;
    reliableMoves(~isnan(reliableMoves)) = cellfun(@(x) max(diff([x(1);x(:)]))<50, preThreshOnsets(~isnan(reliableMoves))) &  ...
        cellfun(@(x) length(unique(x))==1, preThreshDirs(~isnan(reliableMoves)));    
    
    %SANITY CHECK
    responseMade = double(block.events.correctResponseValues(1:length(block.events.feedbackValues))).*(block.events.feedbackValues ~= 0);
    responseMade(block.events.feedbackValues<0) = -1*(responseMade(block.events.feedbackValues<0));
    responseMade = responseMade(block.events.feedbackValues ~= 0)';
    responseMade = ((responseMade>0)+1).*(responseMade~=0);
    if mean(firstMoveDir(reliableMoves==1) == responseMade(reliableMoves==1)) < 0.65 && sum(reliableMoves==1) > 25
        warning('Why are most of the movements not in the same direction as the response?!?')
        keyboard
    end
    aligned.firstMoveTimesDirReliable = [firstMoveIdx./sR firstMoveDir reliableMoves];
    %%
end

if any(contains(fineTune, 'wheelTraceTimeValue'))
    if ~exist('timelinehWeelPosition', 'var')
        timelinehWeelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
        timelinehWeelPosition(timelinehWeelPosition > 2^31) = timelinehWeelPosition(timelinehWeelPosition > 2^31) - 2^32;
    end
    changePoints = strfind(diff([0,timelinehWeelPosition'])==0, [1 0]);
    trialStEnIdx = (trialStEnTimes*sR);
    points2Keep = sort([1 changePoints changePoints+1 length(timelinehWeelPosition) ceil(trialStEnIdx(:,1))'+1, floor(trialStEnIdx(:,2))'-1]);
    aligned.wheelTraceTimeValue = [timelineTime(points2Keep)' timelinehWeelPosition(points2Keep)];
end

rawFields = fields(aligned);
for i = 1:length(rawFields)
    currField = rawFields{i};
    currData = aligned.(currField);
    aligned.(currField) = prc.indexByTrial(trialStEnTimes, currData(:,1), currData);
    emptyIdx = cellfun(@isempty, aligned.(currField));

    if any(strcmp(currField, {'audStimOnOff'; 'visStimOnOff'; 'rewardTimes';'wheelTraceTimeValue'}))
        nColumns = max(cellfun(@(x) size(x,2), aligned.(currField)));
        aligned.(currField)(emptyIdx) = {nan*ones(1,nColumns)};
        aligned.(currField) = cellfun(@single,aligned.(currField), 'uni', 0);
    end
    if any(strcmp(currField, {'audStimPeriodOnOff'; 'visStimPeriodOnOff'; 'firstMoveTimesDirReliable'}))
        nColumns = max(cellfun(@(x) size(x,2), aligned.(currField)));
        aligned.(currField)(emptyIdx) = {nan*ones(1, nColumns)};
        aligned.(currField) = single(cell2mat(aligned.(currField)));
    end
end
if isfield(aligned, 'firstMoveTimesDirReliable')
    aligned.firstMoveTimes = aligned.firstMoveTimesDirReliable(:,1);
    aligned.firstMoveDirection = aligned.firstMoveTimesDirReliable(:,2);
    aligned.firstMoveReliable = aligned.firstMoveTimesDirReliable(:,3);
    aligned = rmfield(aligned, 'firstMoveTimesDirReliable');
end

requiredFields = {'audStimOnOff'; 'visStimOnOff'; 'audStimPeriodOnOff';'visStimPeriodOnOff';'firstMoveTimes';...
                  'firstMoveDirection'; 'firstMoveReliable'; 'rewardTimes';'wheelTraceTimeValue'};
for i = 1:length(requiredFields)
    if ~isfield(aligned, requiredFields{i}); aligned.(requiredFields{i}) = trialStEnTimes(:,2)*0+nan; end
end

aligned.alignment = block.alignment;
correctedBlock = block;
correctedBlock.fineTuned = 1;
end
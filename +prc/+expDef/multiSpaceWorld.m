function [newBlock, newParams, newRaw] = multiSpaceWorld(x)
%% A helper function for multisensoySpaceWorld experimental definition that produces standardised files with useful structures for further analysis.

% Inputs
% x---------------------Structure form the "convertExpFiles" function that contains all the file information

% Outputs
% "newBlock" is the new, compact, and restructured block file with the following fields:
%.subject--------------Name of the mouse
%.expDate--------------Date that the experiment was recorded
%.expNum-----------Session number for experiment
%.rigName--------------Name of the rig where the experiment took place
%.expType--------------Type of the rig where the experiment took place
%.trialStartEnd--------nx2 vector of trial [start end] times relative to the start of the experiment (s)
%.stimPeriodStart------nx1 vector of stimulus period start times relative to the start of the experiment (s)
%.closedLoopStart------nx1 vector of times when the mouse can begin interacting with the stimulus, relative to the stimulus appearing (s)
%.feedback-------------nx1 logical which defines whether the mouse got a reward (1) or not (0)
%.correctResponse------nx1 vector indicating whether the correct response was to turn the wheel left (1) or right (2)
%.responseMade---------nx1 vector responses made by the mouse, whether left (1) or right (2)
%.repeatNum------------nx1 vector indicating the number of times a trial was repeated
%.reactionTime---------nx1 vector of times from closed loop onset until the animal made a response
%.clickDuration--------nx1 vector of click durations (normally the same on every trial) (Hz)
%.clickRate------------nx1 vector of click rate (normally the same on every trial) (Hz)
%.audAmplitude---------nx1 vector of auditory amplitude (does not vary within a trial) (arbitrary units)
%.visContrast----------nx1 vector of visual contrast (does not vary within a trial), as a fraction of 100% contrast.
%.audDiff--------------nx1 vector of the difference in auditory amplitude between left and right (right - left)
%.visDiff--------------nx1 vector of the difference in visual contrast between left and right (right - left)
%.audInitialAzimuth----nx1 vector of the initial azimuth of the auditory stimulus (deg)
%.audInitialAzimuth----nx1 vector of the initial azimuth of the visual stimulus (deg)
%.visAltitude----------nx1 vector of the altitude of the visual contrast stimulus (normally the same on every trial) (deg)
%.visSigma-------------nx2 vector of the sigma values for the visual gabors (normally the same on every trial) (deg)
%.conditions-----------nx1 vector indicating the exact condition index for each trial
%.trialType------------nx1 vector indicating whether a trial was blank (0), aud (1), vis (0), coherent (3), or conflict (4)
%.uniqueConditions-----nxm matrix of unique conditions in the session with columns [audAmplitude visContrast audInitialAzimuth visInitialAzimuth]
%.uniqueConditionsIdx--nxm matrix of the index for each uniqueCondition (used in the "conditions" field)


% ADDED only if there was inactivation during the session
%.galvoType------------nx1 vector indicating galvo moveement, which can be stationary (1) or moving between two sites (2) while the laser is on
%.galvoPosition--------nx2 vector of galvo coordinates used on each trial ([LM axis, AP axis]
%.laserType------------nx1 vector indicating laser type, which can be off (0), unilateral (1), or bilateral (2)
%.laserPower-----------nx1 vector of laser power (mW in more recent experiments)
%.laserSession---------nx1 vector indicating whether had no inactivation (0), unilateral (1), bilateral (2), or both (3)
%.laserOnOff-----------nx2 vector of [on off] times for the laser, relative to trial start.

% "newParams" is the new, compact, and restructured parameters file. Each field will be either one number (same value for for every trial),
% or will have columns equal to the columns in the "numRepeats" field (which is for variables that change with trial type). Fields are:
%.subject----------------------Name of the mouse
%.expDate----------------------Date that the experiment was recorded
%.expNum-----------------------Session number for experiment
%.rigName----------------------Name of the rig where the experiment took place
%.expType----------------------Type of the rig where the experiment took place
%.minutesOnRig-----------------Number of minutes spent on the rig
%.wheelGain--------------------The gain of the wheel (somewhat arbitrary without knowing more about the rotary encoder)
%.numRepeats-------------------Relative frequency of each trial type.
%.galvoType--------------------The galvo type (usually only one value)
%.laserPower-------------------Laser power (mW in later files)
%.laserTypeProportions---------Relative proportion of each laser types 0 vs 1 vs 2 (see above)
%.backgroundNoiseAmplitude-----Amplitude of constant background white noise (arbitraty units)
%.maxRepeatIncorrect-----------Maximum number of times an incorrect trial will be repeated if incorrect.
%.visContrast------------------Contrast of visual stimulus (fraction of 100%)
%.audAmplitude-----------------Amplitude of auditory stimulus (arbitrary units)
%.clickDuration----------------The duration of the stimulus click (s)
%.clickRate--------------------The rate of the stimulus click (Hz)
%.visAltitude------------------Altitude of the visual stimulus (deg relative to zero)
%.visSigma---------------------Sigma values for the visual gabor
%.audInitialAzimuth------------Initial azimuth of the auditory stimulus (deg)
%.audInitialAzimuth------------Initial azimuth of the visual stimulus (deg)
%.openLoopDuration-------------Duration of the open loop phase (when turning the wheel will not move the stimulus)
%.delayAfterIncorrect----------How long to wait after an incorrect trial before starting the next trial (s)
%.laserDuration----------------Duraction of laser stimulation
%.closedLoopOnsetToneAmplitude-Amplitude of closed loop onset tone (arbitraty units)
%.delayAfterCorrect------------How long to wait after an correct trial before starting the next trial (s)
%.rewardSize-------------------Reward size on correct responses (ul)
%.noiseBurstAmplitude----------Amplitude of noise burst on incorrect trials (arbitrary units)
%.noiseBurstDuration-----------Duration of noise burst on incorrect trials (s)
%.stimDuration-----------------Duration of stimulus ("inf" means stim continues until response is made) (s)
%.preStimQuiescentRange--------TODO--define this
%.preStimQuiescentThreshold----TODO--define this
%.rewardTotal------------------Total reward received (ul)
%.totalTrials------------------Total trials performed
%.galvoCoords------------------Set of galvo coordinates used for inactivation
%.numberConditions-------------Number of unique conditions presented
%.audPerformance---------------% correct on auditory trials
%.visPerformance---------------% correct on visual trials
%.mulPerformance---------------% correct on coherent multisensory trials
%.validTrials------------------Number of valid trials (after repeats etc. have been removed)

% "newRaw" is a structure comprising potentially useful raw data (such as wheel movement and timeline data) which is not used for a lot of analyses and
% so should only be loaded if necessary (as it is large). Fields are:
%.visStimOnOff---------{nx1} cell, each contatining an [nx2] vector of [time(s), value(on/off is 1/0)] for the visual stimulus
%.audStimOnOff---------{nx1} cell, each contatining an [nx2] vector of [time(s), value(on/off is 1/0)] for the auditory stimulus
%.visAzimuthTimeValue--{nx1} cell, each contatining an [nx2] vector of [time(s), value(deg)] for the visual stimulus azimth
%.audAzimuthTimeValue--{nx1} cell, each contatining an [nx2] vector of [time(s), value(deg)] for the auditory stimulus azimth
%.rawWheelTimeValue----[nx2] vector of wheel [times(s), position(deg)] over the course of the entire session

%% Convert to shorter names for ease of use later
v = x.standardizedBlock.paramsValues;  %Parameter values at start of trial
e = x.standardizedBlock.events;        %Event structure
n = x.newBlock;                        %newBlock, already populated with subject, expDate, expNum, rigName, and expType
p = x.standardizedParams;              %Parameter values at start of entire session (includes multiple values for different conditions
vIdx = x.validTrials;                  %Indices of valid trials (0 for repeats)

%% Remove excess trials if there are more than 100 total trials (in this case, the mouse was likely still learning)
%Note: we cannot perform this stage before the stage above as it will mess with the calculation of totalRepeats
if sum(vIdx) > 150
    %We remove the first 10 and last 10 correct trials for each session we use -1 because we only want to remove these extra trials from totalRepeats.
    vIdx = double(vIdx);
    stimStartTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1);
    quickResponses = (e.feedbackTimes(1:length(vIdx)) - stimStartTimes(1:length(vIdx)))<1.5;
    vIdx(1:max(find(vIdx==1 & e.responseTypeValues(1:length(vIdx))~=0 & quickResponses, 5, 'first'))) = -1;
    vIdx(min(find(vIdx==1 & e.responseTypeValues(1:length(vIdx))~=0 & quickResponses, 5, 'last')):end) = -1;
    
    %Remove trials in which the laser "trasitionTimes" are more than 90% of the time until the TTL pulse that activates the laser. Otherwise, cannot
    %be confident that the laser was ready to receive the pulse.#
    if p.laserSession
    trasitionTimes = e.laserInitialisationTimes(1:length(vIdx))-e.newTrialTimes(1:length(vIdx));
    timeToTTL = (e.galvoTTLTimes(1:length(vIdx))-e.newTrialTimes(1:length(vIdx)))*0.9;
    vIdx(trasitionTimes(:)>timeToTTL(:) & vIdx(:)==1)=-1;
    
%     if  p.laserPower>0 && p.laserOnsetDelays(1) ~= p.laserOnsetDelays(2)
%         x.timeline = load(x.rawTimeline); x.timeline=x.timeline.Timeline;
%         timelineInfo = prc.extractTimelineInfo(x);
%         vIdx(timelineInfo.badTrials) = -1;
%     end
    end
    vIdx = vIdx>0;
end

%% The number of repeats and timeouts for each trial type presented
%maxRepeatIdx is the set of indices when repeat numbers decrease (i.e. when a maxRepeat is reached, or a repeated trial is performed correctly.
%potentialRepeats is the set of indices for when a trial is incorrect, and it had the potential to repeat.
%totalRepeats is the total number of times each trial was repeated (we subtract 1 so it will be zero if a trial was correct the first time)
timeOutsBeforeResponse = 0*vIdx';
sequentialTimeOuts = 0*vIdx';
sequentialTimeOuts(vIdx &  e.responseTypeValues(1:length(vIdx)) == 0) = 1;
for i = find(vIdx)
    if i < length(vIdx) && e.responseTypeValues(i) == 0
        nextResponse = min([i+find(e.responseTypeValues(i+1:length(vIdx))~=0 | e.repeatNumValues(i+1:length(vIdx))==1,1), length(vIdx)+1]);
        sequentialTimeOuts(i) = nextResponse - i;
        if e.repeatNumValues(nextResponse)==1 || nextResponse >= length(vIdx); continue; end
        vIdx(nextResponse) = 1;
        timeOutsBeforeResponse(nextResponse) = nextResponse-i;
    end
end

repeatsAfterResponse = 0*vIdx';
for i = find(vIdx)
    if i < length(vIdx) && e.responseTypeValues(i) ~= 0
        nextNewTrial = min([i+find(e.repeatNumValues(i+1:length(vIdx))==1,1) length(vIdx)+1]);
        repeatsAfterResponse(i) = nextNewTrial-i-1;
    end
end
%% Populate fields of "n" with basic trial data
%stimPeriodStart extracts times when stimPeriodOnOffValues is 1 (which is when this period starts).
%We also remove times when the first newTrialTime is greater than that the first stimPeriodStart, an error that can occur on the first trial.
%responseTime are the times taken between the stimulus starting and the response being made (including open loop period). Must use
%"1:length(stimPeriodStart)" because if a trial is interupped there can be more stimPeriodStart values than feedback values.
stimPeriodStart = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues == 1)';
stimPeriodStart = stimPeriodStart(vIdx);
feedbackTimes = e.feedbackTimes(vIdx)';
feedbackValues = e.feedbackValues(vIdx)';
timeOuts = feedbackValues==0;
responseTime = feedbackTimes-stimPeriodStart;

if length(p.audAmplitude)==1; p.audAmplitude = repmat(p.audAmplitude,1,length(p.numRepeats)); end    %Make sure there is a value for each condition
if length(p.visContrast)==1; p.visContrast = repmat(p.visContrast,1,length(p.numRepeats)); end       %Make sure there is a value for each condition
audAmplitude = [v(vIdx).audAmplitude]';               %Convert amplitudes to matrix. Assumes one value for each trial.
visContrast = [v(vIdx).visContrast]';                 %Convert amplitudes to matrix. Assumes one value for each trial.
correctResponse = [v(vIdx).correctResponse]';         %Convert correctResponse on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth = [v(vIdx).audInitialAzimuth]';     %Convert audInitialAzimuth on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth(audAmplitude==0) = inf;             %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
p.audInitialAzimuth(p.audAmplitude == 0) = inf;       %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
visInitialAzimuth = [v(vIdx).visInitialAzimuth]';     %Convert visInitialAzimuth on each trial to matrix. Assumes one value for each trial.
visInitialAzimuth(visContrast==0) = inf;              %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)
p.visInitialAzimuth(p.visContrast == 0) = inf;        %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)

%Get trial start/end times for valid trials
trialStartTimes = e.newTrialTimes(vIdx)';
trialEndTimes = e.endTrialTimes(vIdx)';
trialTimes = [trialStartTimes trialEndTimes];

%Process the raw data to be stored separately (because it is large). These are the times at which the visual and auditory stimuli turned on/off and
%all the wheel values. All are indexed by trial using the "indexByTrial" function, and times are relative to stimulus onset.
r.visStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues'], stimPeriodStart, [1 0]);
r.audStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues'], stimPeriodStart, [1 0]);
r.wheelTimeValue = prc.indexByTrial(trialTimes, n.rawWheelTimeValue(:,1), [n.rawWheelTimeValue(:,1) n.rawWheelTimeValue(:,2)], stimPeriodStart, [1 0]);
r.wheelTimeValue(cellfun(@isempty, r.wheelTimeValue)) = deal({[0 0]});
r.wheelTimeValue = cellfun(@(x) [x(:,1) x(:,2)-x(find([x(1:end-1,1);1]>0,1),2)], r.wheelTimeValue, 'uni', 0);
n = rmfield(n, 'rawWheelTimeValue');

%As above, but for the auditory and visual azimuth. These are smaller in size, so we save in the main structure.
r.visAzimuthTimeValue = prc.indexByTrial(trialTimes, e.visAzimuthTimes', [e.visAzimuthTimes' e.visAzimuthValues'], stimPeriodStart, [1 0]);
r.audAzimuthTimeValue = prc.indexByTrial(trialTimes, e.audAzimuthTimes', [e.audAzimuthTimes' e.audAzimuthValues'], stimPeriodStart, [1 0]);

%Get closed loop start times, relative to the stimulus start times (likely to all be the same for a constant delay)
closedLoopStart = e.closedLoopOnOffTimes(e.closedLoopOnOffValues == 1)';
closedLoopStart = closedLoopStart(vIdx) - stimPeriodStart;

%Calculate an approximate time to the first wheel movement. This is different from the response time in that it is based on wheel movement, rather
%than the time when the threshold was reached. WheelMove is an interpolation of the wheel movement (to get it's nearest position at every ms).
%wheelToThresh uses the difference between the wheel position at the closed loop start and threshold to calculate the change in wheel value that
%represents a response (this can be very different for different rotary encoders). timeToWheelMove is then the time at which wheelMove exceeds 25% of
%wheelToThresh
wheelMove = cellfun(@(x) [(-1:0.001:x(end,1))', (interp1(x(diff(x(:,1))~=0,1), x(diff(x(:,1))~=0,2), -1:0.001:x(end,1), 'nearest'))'], r.wheelTimeValue(~timeOuts), 'uni', 0);
wheelToThresh = arrayfun(@(x,y,z) x{1}(x{1}(:,1)>y & x{1}(:,1)<z,2), wheelMove, closedLoopStart(~timeOuts), feedbackTimes(~timeOuts)-stimPeriodStart(~timeOuts), 'uni', 0);
wheelToThresh = nanmedian(abs(cellfun(@(x) x(end)-x(1), wheelToThresh)));
timeToWheelMove = feedbackValues;
timeToWheelMove(~timeOuts) = cell2mat(cellfun(@(x) x(find(abs(x(:,2))>wheelToThresh/4 & x(:,1)>0,1),1), wheelMove, 'uni', 0));

%Get the response the mouse made on each trial based on the correct response and then taking the opposite for incorrect trials. NOTE: this will not
%work for a task with more than two response options.
responseMade = double(correctResponse).*~timeOuts;
responseMade(feedbackValues<0) = -1*(responseMade(feedbackValues<0));

%allConditions is all the conditions the mouse actually performed, where conditions can be completely defined by audAmplitude, visContrast, and the
%intial azimuth of aud and vis stimuli.
%Unique conditions is based on the parameter set, and includes all possible conditions that the mouse could have experienced. We repeat audAmplitude
%because in some cases it will only have one value, and in some cases it will have more;
allConditions = [audAmplitude visContrast audInitialAzimuth visInitialAzimuth];
uniqueConditions = unique([p.audAmplitude' p.visContrast' p.audInitialAzimuth' p.visInitialAzimuth'], 'rows');

%Create a set of unique conditions, where each row is a condition in the order: [zero conditions; right conditions; left conditions].
leftInitialConditions = uniqueConditions(uniqueConditions(:,end-1)< 0 | ((isinf(uniqueConditions(:,end-1)) | ~(uniqueConditions(:,end-1))) & uniqueConditions(:,end)<0),:);
if size(leftInitialConditions,1)~=floor(size(uniqueConditions,1)/2); warning('Why are half conditions not left conditions?'); end
zeroConditions = uniqueConditions(all([any(uniqueConditions(:,[1,3])==0,2) any(uniqueConditions(:,[2,4])==0,2)],2),:);
rightInitialConditions = [leftInitialConditions(:,1:2) leftInitialConditions(:,end-1:end)*-1];
rightInitialConditions(isinf(rightInitialConditions)) = inf;
rightInitialConditions = [rightInitialConditions; uniqueConditions(~ismember(uniqueConditions, [zeroConditions; leftInitialConditions; rightInitialConditions], 'rows'),:)];
uniqueConditions = [zeroConditions; rightInitialConditions; leftInitialConditions];

%For each trial, find which row of leftInitialConditions or rightInitialConditions it belongs to, and check that no rows belong to both. We generate a
%condition index where each trial is positive if right, -ve if left, and the opposite of any conditionIdx is the inverse sign of that conditionIdx. We
%also create uniqueConditionReference which has the corresponding condition for each row of uniqueConditions. Finally, we create conditionRowIdx which (for
%every trial) inditcates which row of the uniqueConditions table that trial corresponds to.
[~, rightConditionsIdx] = ismember(allConditions, rightInitialConditions, 'rows');
[~, leftConditionsIdx] = ismember(allConditions, leftInitialConditions, 'rows');
if any(all([rightConditionsIdx~=0, leftConditionsIdx~=0],2)); error('Detect same condition as being Left and Right'); end
conditionLabel = rightConditionsIdx + -1*leftConditionsIdx;
uniqueConditionRowLabels = [0*find(zeroConditions,1); (1:size(rightInitialConditions,1))'; -1*(1:size(leftInitialConditions,1))'];
[~, conditionRowIdx] = ismember(conditionLabel, uniqueConditionRowLabels);
uniqueDiff = [uniqueConditions(:,3) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];


%Create a "trialType" field which is 0,1,2,3,4 for blank, auditory, visual, coherent, and incoherent trials.
audTrial = (visContrast==0 | visInitialAzimuth==0) & (audAmplitude>0 & audInitialAzimuth~=0);
visTrial = (audAmplitude==0 | audInitialAzimuth==0) & (visContrast>0 & visInitialAzimuth~=0);
coherentTrial = sign(visInitialAzimuth.*audInitialAzimuth)>0 & audAmplitude>0 & visContrast>0;
coflictTrial = sign(visInitialAzimuth.*audInitialAzimuth)<0 & audAmplitude>0 & visContrast>0;
trialType = (zeros(length(coflictTrial),1)+(audTrial+visTrial*2+coherentTrial*3+coflictTrial*4));

audDiff = uniqueDiff(conditionRowIdx, 1);
visDiff = uniqueDiff(conditionRowIdx, 2);
%Create vectors that indicate the separated values for contrast and audio azimuth on left and right (used for modeling)

%Populate n with all fields;
n.trialStartEnd = [trialStartTimes trialEndTimes];
n.stimPeriodStart = stimPeriodStart;
n.closedLoopStart = closedLoopStart;
n.correctResponse = (correctResponse>0)+1;
n.feedback = feedbackValues;
n.responseTime = responseTime;
n.timeToWheelMove = timeToWheelMove;
n.responseMade = ((responseMade>0)+1).*(responseMade~=0);
n.trialType = trialType;
n.sequentialTimeOuts = sequentialTimeOuts(vIdx);
n.timeOutsBeforeResponse = timeOutsBeforeResponse(vIdx);
n.repeatsAfterResponse = repeatsAfterResponse(vIdx);
n.audAmplitude = audAmplitude;
n.audInitialAzimuth = audInitialAzimuth;
n.audDiff = audDiff;
n.audValues = unique(uniqueDiff(:,1));
n.visContrast = visContrast;
n.visInitialAzimuth = visInitialAzimuth;
n.visDiff = visDiff;
n.visValues = unique(uniqueDiff(:,2));
n.uniqueConditions = uniqueConditions;
n.uniqueDiff = uniqueDiff;
n.uniqueConditionRowLabels = uniqueConditionRowLabels;
n.conditionLabelRow = [conditionLabel conditionRowIdx]; 

[n.laserType, n.laserPower, n.galvoType, n.laserOnsetDelay] = deal(n.feedback*nan);
[n.galvoPosition,  n.laserOnOff] = deal(n.feedback*[nan nan]);
p.galvoCoords = [nan nan];
if p.laserSession
    %Galvo position is the position of the galvos on each trial. It is changed so that for bilateral trials, the ML axis is always positive (bilateral
    %trials are when the laserTypeValue for that trial was 2). Note that the galvoPosValues output from the expDef are indices for the galvoCoords (with a
    %-ve index indicating the left himisphere). That's why we need to get the galvo posiiton on each trial by using the abs of this index and then
    %multiplying the ML coordinate by the sign of the original index.
    n.laserType = e.laserTypeValues(vIdx)';
    n.laserPower = (e.laserPowerValues(vIdx)');
    n.galvoType = e.galvoTypeValues(vIdx)';
    
    p.galvoCoords = e.galvoCoordsValues(:,1:2);
    galvoPosValues = e.galvoPosValues(vIdx)';
    galvoPosition = p.galvoCoords(abs(galvoPosValues),:);
    galvoPosition(e.laserTypeValues(vIdx)'~=2,1) = galvoPosition(e.laserTypeValues(vIdx)'~=2,1).*sign(galvoPosValues(e.laserTypeValues(vIdx)'~=2));
    n.galvoPosition = galvoPosition;
   
    if exist('timelineInfo', 'var')
        n.laserOnsetDelay = timelineInfo.laserOnsetDelay(vIdx);
        n.laserOnOff = [timelineInfo.laserOnTime(vIdx), timelineInfo.laserOnTime(vIdx)+[v(vIdx).laserDuration]']-trialStartTimes;
    end
end

%%


%Create some useful grids of aud and vis values. Add these to the block.
[n.grids.visValues, n.grids.audValues] = meshgrid(n.visValues, n.audValues);
[~, gridIdx] = ismember(uniqueDiff, [n.grids.audValues(:) n.grids.visValues(:)], 'rows');
n.grids.conditions = nan*ones(length(n.audValues), length(n.visValues));
n.grids.conditions(gridIdx) = uniqueConditionRowLabels;

if p.laserSession; laserOff = n.laserType == 0; else, laserOff = n.feedback*0+1; end
p.maxRepeatIncorrect = max(p.maxRepeatIncorrect);
p.numberConditions = length(unique([audAmplitude, visContrast audInitialAzimuth visInitialAzimuth], 'rows'));
p.audPerformance = round(mean(n.feedback(n.trialType==1 & laserOff & n.responseMade~=0)>0)*100);
p.visPerformance = round(mean(n.feedback(n.trialType==2 & laserOff & n.responseMade~=0)>0)*100);
p.mulPerformance = round(mean(n.feedback(n.trialType==3 & laserOff & n.responseMade~=0)>0)*100);
p.validResponses = sum(n.responseMade~=0);
p.validTrials = sum(vIdx);

x.validTrials = vIdx;
newParams = p;
newBlock = n;
p.orignal = x.oldParams;

%Remove fields of the raw data where the trial lasts longer than 5s. The data from there trials aren't likely to be useful (since the mouse stayed
%still for a long time) and can take up a lot of space in the mat file.
for i = fields(r)'; newRaw.(i{1}) = r.(i{1}); newRaw.(i{1})(newBlock.responseTime>5) = {[]}; end

%% Check that all expected fields exist
expectedBlockFields ={'audAmplitude';'audDiff';'audInitialAzimuth';'audValues';'closedLoopStart';'conditionLabelRow';'correctResponse';'expDate';...
    'expNum';'expType';'feedback';'grids';'repeatsAfterResponse';'responseMade';'responseTime';'rigName';'sequentialTimeOuts';'stimPeriodStart';...
    'subject';'timeOutsBeforeResponse';'timeToWheelMove';'trialStartEnd';'trialType';'uniqueConditionRowLabels';'uniqueConditions';'uniqueDiff';...
    'visContrast';'visDiff';'visInitialAzimuth';'visValues';'laserType';'laserPower';'galvoType';'galvoPosition';'laserOnOff';'laserOnsetDelay'};

expectedParamFields ={'audAmplitude';'audInitialAzimuth';'audPerformance';'backgroundNoiseAmplitude';'clickDuration';'clickRate';...
    'closedLoopOnsetToneAmplitude';'delayAfterCorrect';'delayAfterIncorrect';'expDate';'expNum';'expType';'galvoCoords';...
    'galvoType';'laserDuration';'laserOnsetDelays';'laserPower';'laserSession';'laserTypeProportions';'maxRepeatIncorrect';'minutesOnRig';...
    'mulPerformance';'noiseBurstAmplitude';'noiseBurstDuration';'numRepeats';'numberConditions';'openLoopDuration';'postQuiescentDelay';...
    'preStimQuiescentRange';'preStimQuiescentThreshold';'responseWindow';'rewardSize';'rewardTotal';'rigName';'stimDuration';'subject';...
    'totalTrials';'validResponses';'validTrials';'visAltitude';'visContrast';'visInitialAzimuth';'visPerformance';'visSigma';'waveformType';'wheelGain'};

if any(~contains(fields(newBlock), expectedBlockFields))
    keybaord;
elseif any(~contains(expectedBlockFields, fields(newBlock)))
    keyboard;
    excessfields = expectedBlockFields(~contains(expectedBlockFields, fields(newBlock)));
    fprintf('Field mistmatch in block file %s - \n', excessfields{:});
end

if any(~contains(fields(newParams), expectedParamFields)) || any(~contains(expectedParamFields, fields(newParams)))
    error('Field mistmatch in param file');
end
end
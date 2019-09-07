function x = multiSpaceWorld(x)
%% A helper function for multisensoySpaceWorld experimental definition that produces standardised files with useful structures for further analysis.

% Inputs
% x---------------------Structure form the "convertExpFiles" function that contains all the file information

% Outputs
% "x.newBlock" is the new, compact, and restructured block file with the following fields:
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


%.galvoType------------nx1 vector indicating galvo moveement, which can be stationary (1) or moving between two sites (2) while the laser is on
%.galvoPosition--------nx2 vector of galvo coordinates used on each trial ([LM axis, AP axis]
%.laserType------------nx1 vector indicating laser type, which can be off (0), unilateral (1), or bilateral (2)
%.laserPower-----------nx1 vector of laser power (mW in more recent experiments)
%.laserSession---------nx1 vector indicating whether had no inactivation (0), unilateral (1), bilateral (2), or both (3)
%.laserOnOff-----------nx2 vector of [on off] times for the laser, relative to trial start.


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
n = x.newBlock;                        %x.newBlock, already populated with subject, expDate, expNum, rigName, and expType
p = x.standardizedParams;              %Parameter values at start of entire session (includes multiple values for different conditions
vIdx = e.repeatNumValues(1:length(x.standardizedBlock.events.endTrialTimes))==1;            %Indices of valid trials (0 for repeats)

%% Trim trials if there are more than 100 total trials (if not, the mouse was likely still learning)
if sum(vIdx) > 150
    %We remove the first 10 and last 10 correct trials for each session we use -1 because we only want to remove these extra trials from totalRepeats.
    vIdx = double(vIdx);
    stimStartTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1);
    quickResponses = (e.feedbackTimes(1:length(vIdx)) - stimStartTimes(1:length(vIdx)))<1.5;
    vIdx(1:max(find(vIdx==1 & e.responseTypeValues(1:length(vIdx))~=0 & quickResponses, 5, 'first'))) = -1;
    vIdx(min(find(vIdx==1 & e.responseTypeValues(1:length(vIdx))~=0 & quickResponses, 5, 'last')):end) = -1;
    
    %Remove trials in which the laser "trasitionTimes" are more than 90% of the time until the TTL pulse that activates the laser. Otherwise, cannot
    %be confident that the laser was ready to receive the pulse.
    if strcmp(n.expType, 'inactivation')
        trasitionTimes = e.laserInitialisationTimes(1:length(vIdx))-e.newTrialTimes(1:length(vIdx));
        timeToTTL = (e.galvoTTLTimes(1:length(vIdx))-e.newTrialTimes(1:length(vIdx)))*0.9;
        vIdx(trasitionTimes(:)>timeToTTL(:) & vIdx(:)==1)=-1;
    end
    vIdx = vIdx>0;
end

%% The number of repeats and timeouts for each trial type presented
%maxRepeatIdx is the set of indices when repeat numbers decrease (i.e. when a maxRepeat is reached, or a repeated trial is performed correctly.
%potentialRepeats is the set of indices for when a trial is incorrect, and it had the potential to repeat.
%totalRepeats is the total number of times each trial was repeated (we subtract 1 so it will be zero if a trial was correct the first time)
for i = find(vIdx)
    if i < length(vIdx) && e.responseTypeValues(i) == 0
        nextResponse = min([i+find(e.responseTypeValues(i+1:length(vIdx))~=0 | e.repeatNumValues(i+1:length(vIdx))==1,1), length(vIdx)+1]);
        if e.repeatNumValues(nextResponse)==1 || nextResponse >= length(vIdx); continue; end
        vIdx(nextResponse) = 1;
    end
end

%% Populate fields of "n" with basic trial data
%stimPeriodStart extracts times when stimPeriodOnOffValues is 1 (which is when this period starts).
%timeToFeedback are the times taken between the stimulus starting and the response being made (including open loop period). Must use
%"1:length(stimPeriodStart)" because if a trial is interupped there can be more stimPeriodStart values than feedback values.
eIdx = 1:length(e.endTrialTimes);
vIdx = vIdx(eIdx);

if length(p.audAmplitude)==1; p.audAmplitude = repmat(p.audAmplitude,1,length(p.numRepeats)); end    %Make sure there is a value for each condition
if length(p.visContrast)==1; p.visContrast = repmat(p.visContrast,1,length(p.numRepeats)); end       %Make sure there is a value for each condition
audAmplitude = [v(eIdx).audAmplitude]';               %Convert amplitudes to matrix. Assumes one value for each trial.
visContrast = [v(eIdx).visContrast]';                 %Convert amplitudes to matrix. Assumes one value for each trial.
correctResponse = [v(eIdx).correctResponse]';         %Convert correctResponse on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth = [v(eIdx).audInitialAzimuth]';     %Convert audInitialAzimuth on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth(audAmplitude==0) = inf;             %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
p.audInitialAzimuth(p.audAmplitude == 0) = inf;       %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
visInitialAzimuth = [v(eIdx).visInitialAzimuth]';     %Convert visInitialAzimuth on each trial to matrix. Assumes one value for each trial.
visInitialAzimuth(visContrast==0) = inf;              %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)
p.visInitialAzimuth(p.visContrast == 0) = inf;        %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)

%Get trial start/end times for valid trials
trialTimes = [e.newTrialTimes(eIdx)' e.endTrialTimes(eIdx)'];
stimPeriodStart = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues == 1)'; 
stimPeriodStart = stimPeriodStart(eIdx);
closedLoopStart = e.closedLoopOnOffTimes(e.closedLoopOnOffValues == 1)'; 
closedLoopStart = closedLoopStart(eIdx);
feedbackTimes = e.feedbackTimes(eIdx)';
feedbackValues = e.feedbackValues(eIdx)';
timeOuts = feedbackValues==0;
timeToFeedback = feedbackTimes-stimPeriodStart;

%Process the raw data to be stored separately (because it is large). These are the times at which the visual and auditory stimuli turned on/off and
%all the wheel values. All are indexed by trial using the "indexByTrial" function, and times are relative to stimulus onset.
times2Subtract = [stimPeriodStart stimPeriodStart*0];
r.visStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues'], times2Subtract);
r.audStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues'], times2Subtract);
r.wheelTimeValue = prc.indexByTrial(trialTimes, n.rawWheelTimeValue(:,1), [n.rawWheelTimeValue(:,1) n.rawWheelTimeValue(:,2)], times2Subtract);
r.wheelTimeValue(cellfun(@isempty, r.wheelTimeValue)) = deal({[0 0]});
r.wheelTimeValue = cellfun(@(x) [x(:,1) x(:,2)-x(find([x(1:end-1,1);1]>0,1),2)], r.wheelTimeValue, 'uni', 0);
n = rmfield(n, 'rawWheelTimeValue');

%As above, but for the auditory and visual azimuth.
r.visAzimuthTimeValue = prc.indexByTrial(trialTimes, e.visAzimuthTimes', [e.visAzimuthTimes' e.visAzimuthValues'], times2Subtract);
r.audAzimuthTimeValue = prc.indexByTrial(trialTimes, e.audAzimuthTimes', [e.audAzimuthTimes' e.audAzimuthValues'], times2Subtract);

%Calculate an approximate time to the first wheel movement. This is different from the response time in that it is based on wheel movement, rather
%than the time when the threshold was reached. WheelMove is an interpolation of the wheel movement (to get it's nearest position at every ms).
%wheelToThresh uses the difference between the wheel position at the closed loop start and threshold to calculate the change in wheel value that
%represents a response (this can be very different for different rotary encoders). timeToWheelMove is then the time at which wheelMove exceeds 25% of
%wheelToThresh
time2closedLoopResponded = closedLoopStart(~timeOuts)-stimPeriodStart(~timeOuts);
timeToFeedbackResponded = timeToFeedback(~timeOuts);
wheelMove = cellfun(@(x) [(-1:0.01:x(end,1))', (interp1(x(diff(x(:,1))~=0,1), x(diff(x(:,1))~=0,2), -1:0.01:x(end,1), 'nearest'))'], r.wheelTimeValue(~timeOuts), 'uni', 0);
wheelToThresh = arrayfun(@(x,y,z) x{1}(x{1}(:,1)>y & x{1}(:,1)<z,2), wheelMove, time2closedLoopResponded, timeToFeedbackResponded, 'uni', 0);
wheelToThresh = nanmedian(abs(cellfun(@(x) x(end)-x(1), wheelToThresh(~cellfun(@isempty, wheelToThresh)))));
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
uniqueConditionLabels = [0*find(zeroConditions,1); (1:size(rightInitialConditions,1))'; -1*(1:size(leftInitialConditions,1))'];
uniqueDiff = [uniqueConditions(:,3) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];
[~, conditionRowIdx] = ismember(conditionLabel, uniqueConditionLabels);


%Create a "trialType" field which is 0,1,2,3,4 for blank, auditory, visual, coherent, and incoherent trials.
trialClass.blank = (visContrast==0 | visInitialAzimuth==0) & (audAmplitude==0 | audInitialAzimuth==0);
trialClass.auditory = (visContrast==0 | visInitialAzimuth==0) & (audAmplitude>0 & audInitialAzimuth~=0);
trialClass.visual = (audAmplitude==0 | audInitialAzimuth==0) & (visContrast>0 & visInitialAzimuth~=0);
trialClass.coherent = sign(visInitialAzimuth.*audInitialAzimuth)>0 & audAmplitude>0 & visContrast>0;
trialClass.conflict = sign(visInitialAzimuth.*audInitialAzimuth)<0 & audAmplitude>0 & visContrast>0;

audDiff = uniqueDiff(conditionRowIdx, 1);
visDiff = uniqueDiff(conditionRowIdx, 2);
%Create vectors that indicate the separated values for contrast and audio azimuth on left and right (used for modeling)

[inactivation.laserType, inactivation.laserPower, inactivation.galvoType, inactivation.laserOnsetDelay] = deal(audDiff*nan);
[inactivation.galvoPosition,  inactivation.laserOnOff] = deal(audDiff*[nan nan]);
if strcmp(n.expType, 'inactivation')
    %Galvo position is the position of the galvos on each trial. It is changed so that for bilateral trials, the ML axis is always positive (bilateral
    %trials are when the laserTypeValue for that trial was 2). Note that the galvoPosValues output from the expDef are indices for the galvoCoords (with a
    %-ve index indicating the left himisphere). That's why we need to get the galvo posiiton on each trial by using the abs of this index and then
    %multiplying the ML coordinate by the sign of the original index.
    inactivation.laserType = e.laserTypeValues(eIdx)';
    inactivation.laserPower = (e.laserPowerValues(eIdx)');
    inactivation.galvoType = e.galvoTypeValues(eIdx)';
    
    galvoCoords = e.galvoCoordsValues(:,1:2);
    galvoPosValues = e.galvoPosValues(eIdx)';
    galvoPosition = galvoCoords(abs(galvoPosValues),:);
    galvoPosition(e.laserTypeValues(eIdx)'~=2,1) = galvoPosition(e.laserTypeValues(eIdx)'~=2,1).*sign(galvoPosValues(e.laserTypeValues(eIdx)'~=2));
    inactivation.galvoPosition = galvoPosition;
   
    if exist('timelineInfo', 'var')
        inactivation.laserOnsetDelay = timelineInfo.laserOnsetDelay(eIdx);
        inactivation.laserOnOff = [timelineInfo.laserOnTime(eIdx), timelineInfo.laserOnTime(eIdx)+[v(eIdx).laserDuration]']-trialStartTimes;
    end
end

%%
%Create some useful grids of aud and vis values. Add these to the block.
[grids.visValues, grids.audValues] = meshgrid(unique(uniqueDiff(:,2)), unique(uniqueDiff(:,1)));
[~, gridIdx] = ismember(uniqueDiff, [grids.audValues(:) grids.visValues(:)], 'rows');
grids.conditionLabels = nan*grids.visValues;
grids.conditionLabels(gridIdx) = uniqueConditionLabels;

if strcmp(n.expType, 'inactivation'); laserOff = inactivation.laserType == 0; else, laserOff = feedbackValues*0+1; end
audPerformance = round(mean(feedbackValues(trialClass.auditory & laserOff & responseMade~=0 & vIdx(:))>0)*100);
visPerformance = round(mean(feedbackValues(trialClass.visual & laserOff & responseMade~=0 & vIdx(:))>0)*100);
mulPerformance = round(mean(feedbackValues(trialClass.coherent & laserOff & responseMade~=0 & vIdx(:))>0)*100);

%Populate n with all fields;
n.performanceAVM = [audPerformance visPerformance mulPerformance];
n.conditionParametersAV = uniqueDiff;
n.conditionLabels = uniqueConditionLabels;
n.trialClass = trialClass; 
n.timings.trialStartEnd = trialTimes;
n.timings.stimPeriodStart = stimPeriodStart;
n.timings.closedLoopStart = closedLoopStart;
n.stim.audAmplitude = audAmplitude;
n.stim.audInitialAzimuth = audInitialAzimuth;
n.stim.audDiff = audDiff;
n.stim.visContrast = visContrast;
n.stim.visInitialAzimuth = visInitialAzimuth;
n.stim.visDiff = visDiff;
n.stim.conditionLabel = conditionLabel; 
n.inactivation = inactivation;
n.outcome.validTrial = vIdx(:);
n.outcome.feedbackGiven = feedbackValues;
n.outcome.timeToFeedback = timeToFeedback;
n.outcome.timeToWheelMove = timeToWheelMove;
n.outcome.responseMade = ((responseMade>0)+1).*(responseMade~=0);
n.grids = grids;
n.params = x.oldParams;

x.newBlock = n;
%Remove fields of the raw data where the trial lasts longer than 5s. The data from there trials aren't likely to be useful (since the mouse stayed
%still for a long time) and can take up a lot of space in the mat file.
for i = fields(r)'; newRaw.(i{1}) = r.(i{1}); newRaw.(i{1})(x.newBlock.outcome.timeToFeedback>5) = {[]}; end
x.newRaw = r;
end
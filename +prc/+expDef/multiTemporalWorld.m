function [newBlock, newParams, newRaw] = multiTemporalWorld(x)
%% A helper function for multiTemporalWorld experimental definition that produces standardised files with useful structures for further analysis.

% Inputs
% x---------------------Structure form the "convertExpFiles" function that contains all the file information

% Outputs
% "newBlock" is the new, compact, and restructured block file with the following fields:
%.subject--------------Name of the mouse
%.expDate--------------Date that the experiment was recorded
%.sessionNum-----------Session number for experiment
%.rigName--------------Name of the rig where the experiment took place
%.rigType--------------Type of the rig where the experiment took place
%.trialStart-----------nx1 vector of trial start times relative to the start of the experiment (s)
%.trialEnd-------------nx1 vector of trial end times relative to the start of the experiment (s)
%.stimPeriodStart------nx1 vector of stimulus period start times relative to the start of the experiment (s)
%.closedLoopStart------nx1 vector of times when the mouse can begin interacting with the stimulus, relative to the stimulus appearing (s)
%.feedback-------------nx1 logical which defines whether the mouse got a reward (1) or not (0)
%.correctResponse------nx1 vector indicating whether the correct response was to turn the wheel left (1) or right (2)
%.response-------------nx1 vector responses made by the mouse, whether left (1) or right (2)
%.repeatNum------------nx1 vector indicating the number of times a trial was repeated
%.reactionTime---------nx1 vector of times from closed loop onset until the animal made a response
%.clickDuration--------nx1 vector of click durations (normally the same on every trial) (Hz)
%.clickRate------------nx1 vector of click rate (normally the same on every trial) (Hz)
%.actualClickRate------nx1 vector the actual click rate may be different to the requested click rate if the requested click rate is non-integer (Hz)
%.audAmplitude---------nx1 vector of auditory amplitude (does not vary within a trial) (arbitrary units)
%.visContrast----------nx1 vector of visual contrast (does not vary within a trial), as a fraction of 100% contrast.
%.audInitialAzimuth----nx1 vector of the initial azimuth of the auditory stimulus (deg)
%.visInitialAzimuth----nx1 vector of the initial azimuth of the visual stimulus (deg)
%.visAltitude----------nx1 vector of the altitude of the visual contrast stimulus (normally the same on every trial) (deg)
%.visSigma-------------nx2 vector of the sigma values for the visual gabors (normally the same on every trial) (deg)
%.galvoType------------nx1 vector indicating galvo moveement, which can be stationary (1) or moving between two sites (2) while the laser is on
%.galvoPosition--------nx2 vector of galvo coordinates used on each trial ([LM axis, AP axis]
%.laserType------------nx1 vector indicating laser type, which can be off (0), unilateral (1), or bilateral (2)
%.laserPower-----------nx1 vector of laser power (mW in more recent experiments)
%.laserSession---------nx1 vector indicating whether had no inactivation (0), unilateral (1), bilateral (2), or both (3)
%.laserOnOff-----------nx2 vector of [on off] times for the laser, relative to trial start.
%.conditions-----------nx1 vector indicating the exact condition index for each trial
%.uniqueConditions-----nxm matrix of unique conditions in the session with columns [audAmplitude visContrast audInitialAzimuth visInitialAzimuth]
%.uniqueConditionsIdx--nxm matrix of the index for each uniqueCondition (used in the "conditions" field)
%.requestedCoherence --The coherence value given to genRandClickTimes, which it will have used to generate the auditory click times
%.actualCoherenceRight-The proportion of auditory clicks which coincide with the clicks of the right visual stimulus
%.actualCoherenceLeft--The proportion of auditory clicks which coincide with the clicks of the left visual stimulus
%
%.clickTimesRight --------------The particular set of randomised click times used for the right visual stimulus for that trial
%.clickTimesLeft ---------------The particular set of randomised click times used for the left visual stimulus for that trial
%.clickTimesAud ----------------The particular set of randomised click times used for the Auditory stimulus for that trial

% "newParams" is the new, compact, and restructured parameters file. Each field will be either one number (same value for for every trial),
% or will have columns equal to the columns in the "numRepeats" field (which is for variables that change with trial type). Fields are:
%.subject----------------------Name of the mouse
%.expDate----------------------Date that the experiment was recorded
%.sessionNum-------------------Session number for experiment
%.rigName----------------------Name of the rig where the experiment took place
%.rigType----------------------Type of the rig where the experiment took place
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
%.visInitialAzimuth------------Initial azimuth of the visual stimulus (deg)
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
%.validTrials------------------Number of valid trials (after repeats etc. have been removed)
%.minInterClickGap-------------The minimum delay required between clicks within the same click train
%.minLeftRightGap--------------The minimum delay required between clicks in the left and right click trains
%.coherentPerformance----------%correct for trials with coherence ~=0.5

% "newRaw" is a structure comprising potentially useful raw data (such as wheel movement and timeline data) which is not used for a lot of analyses and
% so should only be loaded if necessary (as it is large). Fields are:
%.visStimOnOffRightTimeValue---------{nx1} cell, each contatining an [nx2] vector of [time(s), value(on/off is 1/0)] for the right visual stimulus
%.visStimOnOffLeftTimeValue---------{nx1} cell, each contatining an [nx2] vector of [time(s), value(on/off is 1/0)] for the left visual stimulus
%.audStimOnOffTImeValue---------{nx1} cell, each contatining an [nx2] vector of [time(s), value(on/off is 1/0)] for the auditory stimulus
%.visAzimuthRightTimeValue--{nx1} cell, each contatining an [nx2] vector of [time(s), value(deg)] for the right visual stimulus azimuth
%.visAzimuthLeftTimeValue--{nx1} cell, each contatining an [nx2] vector of [time(s), value(deg)] for the left visual stimulus azimuth
%.audAzimuthTimeValue--{nx1} cell, each contatining an [nx2] vector of [time(s), value(deg)] for the auditory stimulus azimth
%.rawWheelTimeValue----[nx2] vector of wheel [times(s), position(deg)] over the course of the entire session

%% Convert to shorter names for ease of use later
v = x.standardizedBlock.paramsValues;  %Parameter values at start of trial
e = x.standardizedBlock.events;        %Event structure
n = x.newBlock;                        %newBlock, already populated with subject, expDate, sessionNum, rigName, and rigType
p = x.standardizedParams;              %Parameter values at start of entire session (includes multiple values for different conditions
vIdx = x.validTrials;                  %Indices of valid trials (0 for repeats)

%% The number of repeats for each trial
%maxRepeatIdx is the set of indices when repeat numbers decrease (i.e. when a maxRepeat is reached, or a repeated trial is performed correctly.
%potentialRepeats is the set of indices for when a trial is incorrect, and it had the potential to repeat. 
%totalRepeats is the total number of times each trial was repeated (we subtract 1 so it will be zero if a trial was correct the first time)
maxRepeatIdx = diff([e.repeatNumValues(1:length(vIdx>0))'; 1])<0;
potentialRepeats = ([v(vIdx>0).maxRepeatIncorrect]>0 & e.feedbackValues(vIdx>0)<=0)';
if potentialRepeats(end) == 1 && vIdx(end)==1; potentialRepeats(end) = 0; end
totalRepeats = potentialRepeats;
totalRepeats(potentialRepeats>0) = e.repeatNumValues(maxRepeatIdx)-1;

%% Remove excess trials if there are more than 100 total trials (in this case, the mouse was likely still learning)
%Note: we cannot perform this stage before the stage above as it will mess with the calculation of totalRepeats 
if sum(vIdx) > 100
    %We remove the first 10 and last 10 correct trials for each session we use -1 because we only want to remove these extra trials from totalRepeats.
    vIdx = double(vIdx);
    vIdx(find(vIdx==1, 10, 'first')) = -1;
    vIdx(find(vIdx==1, 10, 'last')) = -1;
    
    %Remove trials in which the laser "trasitionTimes" are more than 90% of the time until the TTL pulse that activates the laser. Otherwise, cannot
    %be confident that the laser was ready to receive the pulse.
    laserDelays = (1:min([numel(e.galvoTTLTimes) numel(e.laserInitialisationTimes)]))';
    trasitionTimes = e.laserInitialisationTimes(laserDelays);
    timeToTTL = (e.galvoTTLTimes(laserDelays)-e.newTrialTimes(laserDelays))*0.9;
    vIdx(trasitionTimes(:)>timeToTTL(:))=-1;
    
    %Remove the newly eliminated trials from totalRepeats and update vIdx
    totalRepeats(vIdx(vIdx~=0)==-1) = [];
    vIdx = vIdx>0;
end

%% Populate fields of "n" with basic trial data
%stimPeriodStart extracts times when stimPeriodOnOffValues is 1 (which is when this period starts). 
%We also remove times when the first newTrialTime is greater than that the first stimPeriodStart, an error that can occur on the first trial.
%responseTime are the times taken between the stimulus starting and the response being made (including open loop period). Must use 
%"1:length(stimPeriodStart)" because if a trial is interupped there can be more stimPeriodStart values than feedback values. 
stimPeriodStart = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues == 1)'; 
vIdx(e.newTrialTimes(1:length(stimPeriodStart))' > stimPeriodStart) = 0;
stimPeriodStart = stimPeriodStart(vIdx);
feedbackTimes = e.feedbackTimes(vIdx)';
feedbackValues = e.feedbackValues(vIdx)';
caughtMovements = feedbackValues==0; %To identify all trials which ended early due to movement in post stimulus quiescence period
responseTime = feedbackTimes-stimPeriodStart;
responseType = e.responseTypeValues(vIdx)';


if length(p.audAmplitude)==1; p.audAmplitude = repmat(p.audAmplitude,1,length(p.numRepeats)); end    %Make sure there is a value for each condition
if length(p.audStimMobile)==1; p.audStimMobile = repmat(p.audStimMobile,1,length(p.numRepeats)); end %Make sure there is a value for each condition, as it determines whether audAzimuth has any meaning
if length(p.visContrast(1,:))==1; visContrastRight = repmat(p.visContrast(1,:),1,length(p.numRepeats)); else visContrastRight = p.visContrast(1,:); end
if length(p.visContrast(2,:))==1; visContrastLeft = repmat(p.visContrast(2,:),1,length(p.numRepeats)); else visContrastLeft = p.visContrast(2,:); end
audAmplitude = [v(vIdx).audAmplitude]';               %Convert amplitudes to matrix. Assumes one value for each trial.
visContrast = cell2mat({v(vIdx).visContrast}');                
p.galvoCoords = e.galvoCoordsValues(:,1:2);           %Add galvoCoords to the parameter list (take first two columns in case it concatenated across trials)
correctResponse = [v(vIdx).correctResponse]';         %Convert correctResponse on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth = [v(vIdx).audInitialAzimuth]';     %Convert audInitialAzimuth on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth(audAmplitude==0) = inf;             %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
p.audInitialAzimuth(p.audAmplitude == 0) = inf;       %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
audStimMobile = [v(vIdx).audStimMobile]';             %Convert whether the aud stimulus was mobile to a matrix.
visInitialAzimuth = [v(vIdx).visInitialAzimuth]';     %Convert visInitialAzimuth on each trial to matrix. Assumes one value for each trial.
visInitialAzimuth(all(visContrast==0,2)) = inf;              %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)
p.visInitialAzimuth(all(p.visContrast==0,2)) = inf;        %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)
requestedCoherence = e.requestedCoherenceValues(vIdx)';
p.requestedCoherence = (round(p.requestedCoherence.*1000))./1000;

%Galvo position is the position of the galvos on each trial. It is changed so that for bilateral trials, the ML axis is always positive (bilateral
%trials are when the laserTypeValue for that trial was 2. Note that the galvoPosValues output from the expDef are indices for the galvoCoords (with a
%-ve index indicating the left himisphere). That's why we need to get the galvo posiiton on each trial by using the abs of this index and then
%multiplying the ML coordinate by the sign of the original index.
laserTypeValues = e.laserTypeValues(vIdx)';
galvoPosValues = e.galvoPosValues(vIdx)';
galvoPosition = p.galvoCoords(abs(galvoPosValues),:);
galvoPosition(laserTypeValues~=2,1) = galvoPosition(laserTypeValues~=2,1).*sign(galvoPosValues(laserTypeValues~=2));

%Get trial start/end times for valid trials
trialStartTimes = e.newTrialTimes(vIdx)';
trialEndTimes = e.endTrialTimes(vIdx)';
trialTimes = [trialStartTimes trialEndTimes];

%Process the raw data to be stored separately (because it is large). These are the times at which the visual and auditory stimuli turned on/off and
%all the wheel values. All are indexed by trial using the "indexByTrial" function, and times are relative to stimulus onset.
r.visStimOnOffRightTimeValue = prc.indexByTrial(trialTimes, e.visStimOnOffRightTimes', [e.visStimOnOffRightTimes' e.visStimOnOffRightValues'], stimPeriodStart, [1 0]);
r.visStimOnOffLeftTimeValue = prc.indexByTrial(trialTimes, e.visStimOnOffLeftTimes', [e.visStimOnOffLeftTimes' e.visStimOnOffLeftValues'], stimPeriodStart, [1 0]);
r.audStimOnOffTimeValue = prc.indexByTrial(trialTimes, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues'], stimPeriodStart, [1 0]);
r.wheelTimeValue = prc.indexByTrial(trialTimes, n.rawWheelTimeValue(:,1), [n.rawWheelTimeValue(:,1) n.rawWheelTimeValue(:,2)], stimPeriodStart, [1 0]);
r.wheelTimeValue = cellfun(@(x) [x(:,1) x(:,2)-x(find(x(:,1)>0,1),2)], r.wheelTimeValue, 'uni', 0);
n = rmfield(n, 'rawWheelTimeValue');

%As above, but for the auditory and visual azimuth. These are smaller in size, so we save in the main structure.
r.visAzimuthRightTimeValue = prc.indexByTrial(trialTimes, e.visAzimuthRightTimes', [e.visAzimuthRightTimes' e.visAzimuthRightValues'], stimPeriodStart, [1 0]);
r.visAzimuthLeftTimeValue = prc.indexByTrial(trialTimes, e.visAzimuthLeftTimes', [e.visAzimuthLeftTimes' e.visAzimuthLeftValues'], stimPeriodStart, [1 0]);
r.audAzimuthTimeValue = prc.indexByTrial(trialTimes, e.audAzimuthTimes', [e.audAzimuthTimes' e.audAzimuthValues'], stimPeriodStart, [1 0]);

%Get closed loop start times, relative to the stimulus start times (likely to all be the same for a constant delay)
closedLoopStart = e.closedLoopOnOffTimes(e.closedLoopOnOffValues == 1)';
closedLoopStart = prc.indexByTrial(trialTimes, closedLoopStart, closedLoopStart, stimPeriodStart, 0);
closedLoopStart(cellfun(@isempty, closedLoopStart)) = {nan};
closedLoopStart = cell2mat(closedLoopStart) - stimPeriodStart;

%Calculate an approximate time to the first wheel movement. This is different from the response time in that it is based on wheel movement, rather
%than the time when the threshold was reached. WheelMove is an interpolation of the wheel movement (to get it's nearest position at every ms).
%wheelToThresh uses the difference between the wheel position at the closed loop start and threshold to calculate the change in wheel value that
%represents a response (this can be very different for different rotary encoders). timeToWheelMove is then the time at which wheelMove exceeds 25% of
%wheelToThresh
wheelMove = cellfun(@(x) [(-1:0.001:x(end,1))', (interp1(x(:,1), x(:,2), -1:0.001:x(end,1), 'nearest'))'], r.wheelTimeValue, 'uni', 0);
wheelToThresh = arrayfun(@(x,y,z) x{1}(x{1}(:,1)>y & x{1}(:,1)<z,2), wheelMove, closedLoopStart, feedbackTimes - stimPeriodStart, 'uni', 0);
wheelToThresh(responseType==0) = {nan};
wheelToThresh = nanmedian(abs(cell2mat(cellfun(@(x) double(x(end)-x(1)), wheelToThresh, 'uni', 0))));
timeToWheelMove = double(cell2mat(cellfun(@(x) x(find(abs(x(:,2))>wheelToThresh/4 & x(:,1)>0,1),1), wheelMove, 'uni', 0)));

%Get the response the mouse made on each trial based on the correct response and then taking the opposite for incorrect trials. NOTE: this will not
%work for a task with more than two response options.
responseMade = double(correctResponse).*~caughtMovements; %So that early ended trials should be displayed as 0. 
responseMade(feedbackValues<0) = -1*(responseMade(feedbackValues<0));

%allConditions is all the conditions the mouse actually performed, where conditions can be completely defined by audAmplitude, visContrast, the
%intial azimuth of aud and vis stimuli, and the coherence
%Unique conditions is based on the parameter set, and includes all possible conditions that the mouse could have experienced. We repeat audAmplitude
%because in some cases it will only have one value, and in some cases it will have more;
allConditions = [audAmplitude audInitialAzimuth audStimMobile visContrast(:,1) visContrast(:,2) requestedCoherence];
uniqueConditions = unique([p.audAmplitude' p.audInitialAzimuth' p.audStimMobile' visContrastRight' visContrastLeft' p.requestedCoherence'], 'rows');

%Create a set of unique conditions, where each row is a condition in the order: [incoherent conditions; right conditions; left conditions]. 
leftCoherentInitialConditions = uniqueConditions(uniqueConditions(:,end)<0.5,:);
rightCoherentInitialConditions = uniqueConditions(uniqueConditions(:,end)>0.5,:);
incoherentInitialConditions = uniqueConditions(uniqueConditions(:,end)==0.5,:);
if size(leftCoherentInitialConditions,1)~= (size(rightCoherentInitialConditions,1)); warning('Imbalanced right and left conditions'); end
uniqueConditions = [incoherentInitialConditions; rightCoherentInitialConditions; leftCoherentInitialConditions];

%For each trial, find which row of leftCoherentInitialConditions or rightCoherentInitialConditions it belongs to, and check that no rows belong to both. We generate a
%condition index where each trial is positive if right, -ve if left, and the opposite of any conditionIdx is the inverse sign of that conditionIdx. We
%also create uniqueConditionReference which has the corresponding condition for each row of uniqueConditions. Finally, we create conditionRowIdx which (for
%every trial) inditcates which row of the uniqueConditions table that trial corresponds to.

[~, rightConditionsIdx] = ismember(allConditions, rightCoherentInitialConditions, 'rows');
[~, leftConditionsIdx] = ismember(allConditions, leftCoherentInitialConditions, 'rows'); %%%Does not seem to be adding values...
if any(all([rightConditionsIdx~=0, leftConditionsIdx~=0],2)); error('Detect same condition as being Left and Right'); end
conditionLabel = rightConditionsIdx + -1*leftConditionsIdx;
uniqueConditionRowLabels = [0*find(incoherentInitialConditions,1); (1:size(rightCoherentInitialConditions,1))'; -1*(1:size(leftCoherentInitialConditions,1))'];
[~, conditionRowIdx] = ismember(conditionLabel, uniqueConditionRowLabels);

%Populate n with all fields;
n.trialStartEnd = [trialStartTimes trialEndTimes];
n.StimPeriodStart = stimPeriodStart;
n.closedLoopStart = closedLoopStart;
n.rewardAvailable = e.rewardAvailableValues(vIdx)'>0;
n.correctResponse = (correctResponse>0)+1;
n.feedback = feedbackValues; %Changed from feedbackValues>0 so that early ended trials will be notable. Need to change analysis code probably though
n.responseTime = responseTime;
n.timeToWheelMove = timeToWheelMove;
n.responseMade = (responseMade>0)+1.*(responseMade~=0); %Added so that when no response has been made the value will be 0.
n.responseType = responseType;
n.audAmplitude = audAmplitude;
n.audInitialAzimuth = audInitialAzimuth;
n.audStimMobile = audStimMobile;
n.visContrast = visContrast; 
n.visInitialAzimuth = visInitialAzimuth;
n.visAltitude = [v(vIdx).visAltitude]';
n.visSigma = [v(vIdx).visSigma]';
n.galvotype = e.galvoTypeValues(vIdx)';
n.galvoPosition = galvoPosition;
n.laserType = e.laserTypeValues(vIdx)';
n.laserPower = (e.laserPowerValues(vIdx)');
n.laserSession = sum(unique(n.laserType))+(n.laserPower*0);
n.laserOnOff = [e.galvoTTLTimes(vIdx)' e.galvoAndLaserEndTimes(vIdx)']-trialStartTimes;
n.rewardAvailable = e.rewardAvailableValues(vIdx)'>0;
n.totalRepeats = totalRepeats;
n.uniqueConditions = uniqueConditions;
n.uniqueConditionRowLabels = uniqueConditionRowLabels;
n.conditionLabel = conditionLabel;
n.conditionRowIdx = conditionRowIdx;
n.requestedCoherence = requestedCoherence; 
n.postStimQuiescenceDuration = e.postStimQuiescentDurationValues(vIdx)';

numOfClicks = round(unique([v(vIdx).clickRate])*(unique([v(vIdx).stimDuration])));
clickTimesRight = reshape(e.clickTimesValues(1,:), [numOfClicks,length(e.trialNumValues)])';
n.clickTimesRight = clickTimesRight(vIdx,:);
clickTimesLeft = reshape(e.clickTimesValues(2,:), [numOfClicks,length(e.trialNumValues)])';
n.clickTimesLeft = clickTimesLeft(vIdx,:);
clickTimesAud = reshape(e.clickTimesValues(3,:), [numOfClicks,length(e.trialNumValues)])';
n.clickTimesAud = clickTimesAud(vIdx,:);
n.actualClickRate = numOfClicks/unique([v(vIdx).stimDuration]);

% Calculate the actual coherence of the click train
actualCoherenceRight = zeros(length(vIdx(vIdx)),1);
actualCoherenceLeft = zeros(length(vIdx(vIdx)),1);
for i = 1:length(vIdx(vIdx))
    actualCoherenceRight(i,1) = mean(ismember(n.clickTimesRight(i,2:end), n.clickTimesAud(i,2:end)));
    actualCoherenceLeft(i,1) = mean(ismember(n.clickTimesLeft(i,2:end), n.clickTimesAud(i,2:end)));
end
n.actualCoherenceRight = actualCoherenceRight;
n.actualCoherenceLeft = actualCoherenceLeft;

p.maxRepeatIncorrect = max(p.maxRepeatIncorrect);
p.numberConditions = length(unique([audAmplitude audInitialAzimuth visContrast(:,1) visContrast(:,2) requestedCoherence], 'rows')); 
p.coherentPerformance = round(mean(n.feedback(n.requestedCoherence ~=0.5 & n.laserType == 0))*100); 
p.validTrials = sum(vIdx);

x.validTrials = vIdx;
newParams = p;
newBlock = n;

%Remove fields of the raw data where the trial lasts longer than 5s. The data from there trials aren't likely to be useful (since the mouse stayed
%still for a long time) and can take up a lot of space in the mat file.
for i = fields(r)'; newRaw.(i{1}) = r.(i{1}); newRaw.(i{1})(newBlock.responseTime>10) = {[]}; end

%%% The sections on trialType, uniqueDiff, and making grids have been removed

%% Check that all expected fields exist
blockFields = {'subject'; 'expDate';'sessionNum';'rigName';'rigType';'trialStartEnd';'StimPeriodStart';'closedLoopStart';'rewardAvailable'; ...
    'correctResponse';'feedback';'responseTime';'timeToWheelMove';'responseMade';'audAmplitude';'audInitialAzimuth';...
    'visContrast';'visInitialAzimuth';'audStimMobile';'visAltitude';'visSigma';'galvotype'; 'responseType';...
    'galvoPosition';'laserType';'laserPower';'laserSession';'laserOnOff';'totalRepeats';'uniqueConditions';'conditionLabel'; ...
    'conditionRowIdx'; 'uniqueConditionRowLabels';'requestedCoherence';'clickTimesRight';'clickTimesLeft';'clickTimesAud';'actualCoherenceRight'; ...
    'actualCoherenceLeft';'actualClickRate';'postStimQuiescenceDuration'};

prmFields =  {'subject';'expDate';'sessionNum';'rigName';'rigType';'wheelGain';'galvoType';'laserPower';'laserTypeProportions';'backgroundNoiseAmplitude';'maxRepeatIncorrect' ...
    ;'visContrast';'audAmplitude';'clickDuration';'clickRate';'visAltitude';'visSigma';'audInitialAzimuth';'visInitialAzimuth';'openLoopDuration' ...
    ;'delayAfterIncorrect';'laserDuration'; 'closedLoopOnsetToneAmplitude';'delayAfterCorrect';'rewardSize';'noiseBurstAmplitude' ...
    ;'noiseBurstDuration';'stimDuration';'preStimQuiescentRange';'preStimQuiescentThreshold';'rewardTotal'; 'responseWindow' ...
    ;'totalTrials';'minutesOnRig';'galvoCoords';'numberConditions';'coherentPerformance';'validTrials' ...
    ; 'numRepeats'; 'audStimMobile'; 'requestedCoherence'; 'minInterClickGap'; 'minLeftRightGap'; 'postStimQuiescentThreshold'; 'postStimQuiescentDuration'};

if any(~contains(fields(newBlock), blockFields)) || any(~contains(blockFields, fields(newBlock)))
    error('Field mistmatch in block file');
end

if any(~contains(fields(newParams), prmFields)) || any(~contains(prmFields, fields(newParams)))
    error('Field mistmatch in param file');
end
end
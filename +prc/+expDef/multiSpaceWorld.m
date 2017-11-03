function [newBlock, newParams, newRaw] = multiSpaceWorld(x)
%% A helper function for multisensoySpaceWorld experimental definition that produces standardised files with useful structures for further analysis.

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
%.audAmplitude---------nx1 vector of auditory amplitude (does not vary within a trial) (arbitrary units)
%.visContrast----------nx1 vector of visual contrast (does not vary within a trial), as a fraction of 100% contrast.
%.audDiff--------------nx1 vector of the difference in auditory amplitude between left and right (right - left)
%.visDiff--------------nx1 vector of the difference in visual contrast between left and right (right - left)
%.audInitialAzimuth----nx1 vector of the initial azimuth of the auditory stimulus (deg)
%.audInitialAzimuth----nx1 vector of the initial azimuth of the visual stimulus (deg)
%.visAltitude----------nx1 vector of the altitude of the visual contrast stimulus (normally the same on every trial) (deg)
%.visSigma-------------nx2 vector of the sigma values for the visual gabors (normally the same on every trial) (deg)
%.galvoType------------nx1 vector indicating galvo moveement, which can be stationary (1) or moving between two sites (2) while the laser is on
%.galvoPosition--------nx2 vector of galvo coordinates used on each trial ([LM axis, AP axis]
%.laserType------------nx1 vector indicating laser type, which can be off (0), unilateral (1), or bilateral (2)
%.laserPower-----------nx1 vector of laser power (mW in more recent experiments)
%.laserSession---------nx1 vector indicating whether had no inactivation (0), unilateral (1), bilateral (2), or both (3)
%.laserOnOff-----------nx2 vector of [on off] times for the laser, relative to trial start.
%.conditions-----------nx1 vector indicating the exact condition index for each trial
%.trialType------------nx1 vector indicating whether a trial was blank (0), aud (1), vis (0), coherent (3), or conflict (4)
%.uniqueConditions-----nxm matrix of unique conditions in the session with columns [audAmplitude visContrast audInitialAzimuth visInitialAzimuth]
%.uniqueConditionsIdx--nxm matrix of the index for each uniqueCondition (used in the "conditions" field)

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
n = x.newBlock;                        %newBlock, already populated with subject, expDate, sessionNum, rigName, and rigType
p = x.standardizedParams;              %Parameter values at start of entire session (includes multiple values for different conditions
vIdx = x.validTrials;                  %Indices of valid trials (0 for repeats)

%% The number of repeats for each trial
%maxRepeatIdx is the set of indices when repeat numbers decrease (i.e. when a maxRepeat is reached, or a repeated trial is performed correctly.
%potentialRepeats is the set of indices for when a trial is incorrect, and it had the potential to repeat. 
%totalRepeats is the total number of times each trial was repeated (we subtract 1 so it will be zero if a trial was correct the first time)
maxRepeatIdx = diff([e.repeatNumValues(1:length(vIdx>0))'; 1])<0;
potentialRepeats = ([v(vIdx>0).maxRepeatIncorrect]>0 & e.feedbackValues(vIdx>0)<0)';
if potentialRepeats(end) == 1 && vIdx(end)==1; potentialRepeats(end) = 0; end
totalRepeats = potentialRepeats;
totalRepeats(potentialRepeats>0) = e.repeatNumValues(maxRepeatIdx)-1;

%% Remove excess trials if there are more than 100 total trials (in this case, the mouse was likely still learning)
%Note: we cannot perfrom this stage befor the stage above as it will mess with the calculation of totalRepeats 
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
    totalRepeats(vIdx(vIdx==1)==-1) = [];
    vIdx = vIdx~=0;
end

%% Populate fields of "n" with basic trial data
%stimPeriodStart extracts times when stimPeriodOnOffValues is 1 (which is when this period starts). 
%We also remove times when the first newTrialTime is greater than that the first stimPeriodStart, an error that can occur on the first trial.
%responseTime are the times taken between the stimulus starting and the response being made (including open loop period). Must use 
%"1:numel(e.endTrialTimes)" because if a trial is interupped there can be more stimPeriodStart values than feedback values. 
stimPeriodStart = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues == 1)'; 
vIdx(e.newTrialTimes(1:length(stimPeriodStart))' > stimPeriodStart) = 0;
stimPeriodStart = stimPeriodStart(vIdx);
feedbackTimes = e.feedbackTimes(vIdx)';
feedbackValues = e.feedbackValues(vIdx)';
responseTime = feedbackTimes-stimPeriodStart;

%Galvo position is the position of the galvos on each tria. It is changed so that for bilateral trials, the ML axis is always positive
laserTypeValues = e.laserTypeValues(vIdx)';
galvoPosValues = galvoPosValues(vIdx)';
galvoPosition = p.galvoCoords(abs(e.galvoPosValues(vIdx))',:);
galvoPosition(:,1) = galvoPosition(:,1).*sign(e.galvoPosValues(vIdx))';
galvoPosition(e.laserTypeValues(vIdx)'==2,1) = abs(galvoPosition(e.laserTypeValues(vIdx)'==2,1));

audAmplitude = [v.audAmplitude]';               %Convert amplitudes to matrix. Assumes one value for each trial.
visContrast = [v.visContrast]';                 %Convert amplitudes to matrix. Assumes one value for each trial.
p.galvoCoords = unique(e.galvoCoordsValues, 'rows');            %Add galvoCoords to the parameter list
correctResponse = [v.correctResponse]';         %Convert correctResponse on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth = [v.audInitialAzimuth]';     %Convert audInitialAzimuth on each trial to matrix. Assumes one value for each trial.
audInitialAzimuth(audAmplitude==0) = inf;       %Change case when audAmplitude was 0 to have infinite azimuth (an indication it has no azimuth value)
p.audInitialAzimuth(p.audAmplitude == 0) = inf; %Change case when audAmplitude was 0 to have infinite azimuth (an indication it has no azimuth value)
visInitialAzimuth = [v.visInitialAzimuth]';     %Convert visInitialAzimuth on each trial to matrix. Assumes one value for each trial.
visInitialAzimuth(visContrast==0) = inf;        %Change case when visContrast was 0 to have infinite azimuth (an indication it has no azimuth value)
p.visInitialAzimuth(p.visContrast == 0) = inf;  %Change case when visContrast was 0 to have infinite azimuth (an indication it has no azimuth value)

%Get trial start/end times for valid trials
trialStartTimes = e.newTrialTimes(vIdx)';
trialEndTimes = e.endTrialTimes(vIdx)';
trialTimes = [trialStartTimes trialEndTimes];

%Process the raw data to be stored separately (because it is large). These are the times at which the visual and auditory stimuli turned on/off and
%all the wheel values. All are inexed by trial using the "indexByTrial" function, and times are relative to stimulus onset.
r.visStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues'], stimPeriodStart, [1 0]);
r.audStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues'], stimPeriodStart, [1 0]);
r.wheelTimeValue = prc.indexByTrial(trialTimes, n.rawWheelTimeValue(:,1), [n.rawWheelTimeValue(:,1) n.rawWheelTimeValue(:,2)], stimPeriodStart, [1 0]);
r.wheelTimeValue = cellfun(@(x) [x(:,1) x(:,2)-x(find(x(:,1)>0,1),2)], r.wheelTimeValue, 'uni', 0);
n = rmfield(n, 'rawWheelTimeValue');

%As above, but for the auditory and visual azimuth. These are smaller in size, so we save in the main structure.
r.visAzimuthTimeValue = prc.indexByTrial(trialTimes, e.visAzimuthTimes', [e.visAzimuthTimes' e.visAzimuthValues'], stimPeriodStart, [1 0]);
r.audAzimuthTimeValue = prc.indexByTrial(trialTimes, e.audAzimuthTimes', [e.audAzimuthTimes' e.audAzimuthValues'], stimPeriodStart, [1 0]);

%Get closed loop start times, relative to the stimulus start times (likely to all be the same for a constant delay)
closedLoopStart = e.closedLoopOnOffTimes(e.closedLoopOnOffValues == 1)';
closedLoopStart = closedLoopStart(vIdx) - stimPeriodStart;

%Calculate an approximate time to the fisrt wheel movement. This is different from the response time in that it is based on wheel movement, rather
%than the time when the threshold was reached. WheelMove is an interpolation of the wheel movement (to get it's nearest position at every ms).
%wheelToThresh uses the difference between the wheel position at the closed loop start and threshold to calculate the change in wheel value that
%represents a response (this can ve different for different rotary encoders). timeToWheelMove is then the time at which wheelMove exceeds 25% of
%wheelToThresh
wheelMove = cellfun(@(x) [(-1:0.001:x(end,1))', (interp1(x(:,1), x(:,2), -1:0.001:x(end,1), 'nearest'))'], r.wheelTimeValue, 'uni', 0);
wheelToThresh = arrayfun(@(x,y,z) x{1}(x{1}(:,1)>y & x{1}(:,1)<z,2), wheelMove, closedLoopStart, feedbackTimes - stimPeriodStart, 'uni', 0);
wheelToThresh = nanmedian(abs(cellfun(@(x) x(end)-x(1), wheelToThresh)));
timeToWheelMove = double(cell2mat(cellfun(@(x) x(find(abs(x(:,2))>wheelToThresh/4 & x(:,1)>0,1),1), wheelMove, 'uni', 0)));

%Get the response the mouse made on each trial based on the correct response and then taking the opposite for incoorect trials. NOTE: this will not
%work for a task with more than two response options.
responseMade = double(correctResponse);
responseMade(feedbackValues<0) = -1*(responseMade(feedbackValues<0)-2)+1;

n.feedback = feedbackValues>0;
n.correctResponse = uint8(correctResponse(vIdx)>0)+1;

n.trialStart = trialStartTimes;
n.trialEnd = trialEndTimes;
n.stimPeriodStart = stimPeriodStart;
n.closedLoopStart = closedLoopStart;
n.response = response;
n.totalRepeats = totalRepeats';
n.rewardAvailable = e.rewardAvailableValues(vIdx)'>0;
n.responseTime = (responseTime(vIdx));
n.timeToWheelMove = timeToWheelMove;
n.clickDuration = ([v(vIdx).clickDuration]');
n.clickRate = ([v(vIdx).clickRate]');
n.audAmplitude = audAmplitude(vIdx,:);
n.visContrast = visContrast(vIdx,:);
n.audDiff = [];
n.visDiff = [];
n.audInitialAzimuth = (audInitialAzimuth(vIdx));
n.visInitialAzimuth = (visInitialAzimuth(vIdx));
n.visAltitude = ([v(vIdx).visAltitude]');
n.visSigma = ([v(vIdx).visSigma]');
n.galvoType = int8(e.galvoTypeValues(vIdx)');
n.galvoPosition = galvoPosition;
n.laserType = int8(e.laserTypeValues(vIdx)');
n.laserPower = (e.laserPowerValues(vIdx)');
n.laserSession = sum(unique(n.laserType))+(n.laserPower*0);
n.laserOnOff = [e.galvoTTLTimes(vIdx)' e.galvoAndLaserEndTimes(vIdx)']-n.trialStart;

allConditions = [n.audAmplitude n.visContrast n.audInitialAzimuth n.visInitialAzimuth];
uniqueConditions = unique([p.audAmplitude' p.visContrast' p.audInitialAzimuth' p.visInitialAzimuth'], 'rows');

if mod(size(allConditions, 2),4) ~=0; error('Must be an initial azimuth for every amplitude and contrast'); end

tempDat = uniqueConditions; tempDat(isinf(uniqueConditions)) = 0;
rightInitialConditions = uniqueConditions(tempDat(:,end-1)> 0 | (tempDat(:,end-1)==0 & tempDat(:,end)>0),:);
leftInitialConditions = [rightInitialConditions(:,1:2) rightInitialConditions(:,end-1:end)*-1];
leftInitialConditions(isinf(leftInitialConditions)) = inf;

[~, rightConditions] = ismember(allConditions, rightInitialConditions, 'rows');
[~, leftConditions] = ismember(allConditions, leftInitialConditions, 'rows');
if any(all([rightConditions~=0, leftConditions~=0],2)); error('Detect same condition as being Left and Right'); end
if size(rightInitialConditions,1)~=size(leftInitialConditions,1); disp('WARNING: Seem to have asymmetical conditions'); end
n.conditions = (rightConditions + -1*leftConditions);
uniqueConditions = uniqueConditions(~ismember(uniqueConditions, [rightInitialConditions; leftInitialConditions], 'rows'),:);
uniqueConditions = [uniqueConditions; rightInitialConditions; leftInitialConditions];

audTrial = (n.visContrast==0 | n.visInitialAzimuth==0) & (n.audAmplitude>0 & n.audInitialAzimuth~=0);
visTrial = (n.audAmplitude==0 | n.audInitialAzimuth==0) & (n.visContrast>0 & n.visInitialAzimuth~=0);
coherentTrial = sign(n.visInitialAzimuth.*n.audInitialAzimuth)>0 & n.audAmplitude>0 & n.visContrast>0;
coflictTrial = sign(n.visInitialAzimuth.*n.audInitialAzimuth)<0 & n.audAmplitude>0 & n.visContrast>0;
blankTrial = (n.audAmplitude == 0 | n.audInitialAzimuth == 0) & (n.visContrast == 0 | n.visInitialAzimuth == 0);
n.trialType = (~blankTrial.*(audTrial+visTrial*2+coherentTrial*3+coflictTrial*4));

n.uniqueConditions = (uniqueConditions);
if length(unique(uniqueConditions(:,1)))==1
    uniqueDiff = [uniqueConditions(:,3) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];
    if any(uniqueConditions(:,1)); n.audType = 'Azimuth'; else; n.audType = 'None'; end
else 
    uniqueDiff = [uniqueConditions(:,1).*sign(uniqueConditions(:,3)) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];
    n.audType = 'Amplitude';
end
n.uniqueDiff = (uniqueDiff);
n.uniqueConditionsIdx = ([0*find(~any(n.uniqueDiff,2)); (1:size(rightInitialConditions,1))';  -1*(1:size(leftInitialConditions,1))']);
[~, n.conditionsIdx] = ismember(n.conditions, n.uniqueConditionsIdx);
n.conditionsIdx = (n.conditionsIdx);

n.visContrastLeftRight = [n.visContrast.*double(n.visInitialAzimuth<0), double(n.visContrast.*double(n.visInitialAzimuth>0))];
n.audAzimuthLeftRight = [abs(n.audInitialAzimuth.*double(n.audInitialAzimuth<0)) abs(n.audInitialAzimuth.*double(n.audInitialAzimuth>0))];

n.audDiff = n.uniqueDiff(n.conditionsIdx, 1);
% add visDiff here!
n.audValues = (unique(uniqueDiff(:,1)));
n.visValues = (unique(uniqueDiff(:,2)));
[n.grids.visValues, n.grids.audValues] = meshgrid(n.visValues, n.audValues);
[~, gridIdx] = ismember(uniqueDiff, [n.grids.audValues(:) n.grids.visValues(:)], 'rows');
n.grids.conditions = nan*ones(length(n.audValues), length(n.visValues));
n.grids.conditions(gridIdx) = n.uniqueConditionsIdx;

p.maxRepeatIncorrect = max(p.maxRepeatIncorrect);
p.numberConditions = length(unique([audAmplitude, visContrast audInitialAzimuth visInitialAzimuth], 'rows'));
p.audPerformance = round(mean(n.feedback(n.trialType==1 & n.laserType == 0))*100);
p.visPerformance = round(mean(n.feedback(n.trialType==2 & n.laserType == 0))*100);
p.mulPerformance = round(mean(n.feedback(n.trialType==3 & n.laserType == 0))*100);
p.validTrials = sum(vIdx);

x.validTrials = vIdx;
newParams = p;
newBlock = n;
for i = fields(r)'; newRaw.(i{1}) = r.(i{1}); newRaw.(i{1})(newBlock.responseTime>5) = {[]}; end

%% Check that all expected fields exist
blockFields = {'subject';'expDate';'sessionNum';'rigName';'rigType';'trialStart';'trialEnd';'stimPeriodStart';'closedLoopStart';'feedback'; 'conditionsIdx';  ...
    'correctResponse';'responseMade';'repeatNum';'timeToWheelMove';'responseTime';'clickDuration';'clickRate';'audAmplitude';'visContrast';'audDiff';'visDiff'; 'visContrastLeftRight';...
    'audInitialAzimuth';'visInitialAzimuth';'visAltitude';'visSigma';'galvoType';'galvoPosition';'laserType';'laserPower';'laserSession'; 'audAzimuthLeftRight';...
    'laserOnOff';'conditions';'trialType';'uniqueConditions'; 'rewardAvailable'; 'visValues'; 'audValues'; 'grids'; 'uniqueDiff'; 'audType'};
prmFields =  {'subject';'expDate';'sessionNum';'rigName';'rigType';'wheelGain';'galvoType';'laserPower';'laserTypeProportions';'backgroundNoiseAmplitude';'maxRepeatIncorrect' ...
    ;'visContrast';'audAmplitude';'clickDuration';'clickRate';'visAltitude';'visSigma';'audInitialAzimuth';'visInitialAzimuth';'openLoopDuration' ...
    ;'delayAfterIncorrect';'laserDuration'; 'closedLoopOnsetToneAmplitude';'delayAfterCorrect';'rewardSize';'noiseBurstAmplitude' ...
    ;'noiseBurstDuration';'stimDuration';'preStimQuiescentRange';'preStimQuiescentThreshold';'rewardTotal' ...
    ;'totalTrials';'minutesOnRig';'galvoCoords';'numberConditions';'audPerformance';'visPerformance';'mulPerformance';'validTrials'; 'numRepeats'};

if any(~contains(fields(newBlock), blockFields)) || any(~contains(blockFields, fields(newBlock)))
    error('Field mistmatch in block file');
end

if any(~contains(fields(newParams), prmFields)) || any(~contains(prmFields, fields(newParams)))
    error('Field mistmatch in param file');
end
end
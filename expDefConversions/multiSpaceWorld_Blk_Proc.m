function [newBlock, newParams, newRaw] = multiSpaceWorld_Blk_Proc(x)
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

%%
v = x.standardizedBlock.paramsValues;
e = x.standardizedBlock.events;
n = x.newBlock;
p = x.standardizedParams;
tIdx = x.validTrials;

stimPeriodStart = x.standardizedBlock.events.stimPeriodOnOffTimes(x.standardizedBlock.events.stimPeriodOnOffValues == 1)';
reactionTime = e.feedbackTimes(1:length(e.endTrialTimes))'-stimPeriodStart(1:length(e.endTrialTimes));

repeatIdx = diff([x.standardizedBlock.events.repeatNumValues(1:length(tIdx>0)) 1])<0;
repeatNum = int8([x.standardizedBlock.paramsValues(tIdx>0).maxRepeatIncorrect]>0 ...
    & x.standardizedBlock.events.feedbackValues(tIdx>0)<0);
if repeatNum(end) == 1 && tIdx(end)==1; repeatNum(end) = 0; end
repeatNum(repeatNum>0) = x.standardizedBlock.events.repeatNumValues(repeatIdx)-1;

if sum(tIdx) > 100
    newVTrials = double(tIdx);
    newVTrials(find(tIdx==1, 10, 'first')) = 2;
    newVTrials(find(tIdx==1, 10, 'last')) = 2;
    
    laserDelays = min([numel(e.galvoTTLTimes) numel(e.laserInitialisationTimes)]);
    newVTrials(e.laserInitialisationTimes(1:laserDelays)>((e.galvoTTLTimes(1:laserDelays)-e.newTrialTimes(1:laserDelays))*0.9)')=2;
    newVTrials = newVTrials.*tIdx;
    
    repeatNum(newVTrials(tIdx==1)==2) = [];
    tIdx = newVTrials==1;
end
n.trialStart = single(x.standardizedBlock.events.newTrialTimes(tIdx)');
n.trialEnd = x.standardizedBlock.events.endTrialTimes(tIdx)';

stimPeriodStart = indexByTrial(n, stimPeriodStart, stimPeriodStart, 0);
n.stimPeriodStart = single(cell2mat(stimPeriodStart));

audAmplitude = single(cell2mat({v.audAmplitude}'));
visContrast = single(cell2mat({v.visContrast}'));

if size(p.audInitialAzimuth,2) < length(p.numRepeats); p.audInitialAzimuth = repmat(p.audInitialAzimuth,1,length(p.numRepeats)); end
if size(p.audAmplitude,2) < length(p.numRepeats); p.audAmplitude = repmat(p.audAmplitude,1,length(p.numRepeats)); end
if size(p.visInitialAzimuth,2) < length(p.numRepeats); p.visInitialAzimuth = repmat(p.visInitialAzimuth,1,length(p.numRepeats)); end
if size(p.visContrast,2) < length(p.numRepeats); p.visContrast = repmat(p.visContrast,1,length(p.numRepeats)); end

p.galvoCoords = e.galvoCoordsValues(:,1:2);
correctResponse = [v.correctResponse]';  %Correct Response
audInitialAzimuth = [v.audInitialAzimuth]'; audInitialAzimuth(audAmplitude==0) = inf;
p.audInitialAzimuth(p.audAmplitude == 0) = inf;
visInitialAzimuth = [v.visInitialAzimuth]'; visInitialAzimuth(visContrast==0) = inf;
p.visInitialAzimuth(p.visContrast == 0) = inf;


audDiff = audAmplitude.*(audInitialAzimuth>0)-audAmplitude.*(audInitialAzimuth<0);
visDiff = visContrast.*(visInitialAzimuth>0)-visContrast.*(visInitialAzimuth<0);

r.visStimOnOff = indexByTrial(n, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues'], [1 0]);
r.audStimOnOff = indexByTrial(n, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues'], [1 0]);

r.visAzimuthTimeValue = indexByTrial(n, e.visAzimuthTimes', [e.visAzimuthTimes' e.visAzimuthValues'], [1 0]);
r.audAzimuthTimeValue = indexByTrial(n, e.audAzimuthTimes', [e.audAzimuthTimes' e.audAzimuthValues'], [1 0]);

closedLoopStart = e.closedLoopOnOffTimes(e.closedLoopOnOffValues == 1)';
closedLoopStart = indexByTrial(n, closedLoopStart, closedLoopStart, 1);
n.closedLoopStart = single(cell2mat(closedLoopStart));

n.feedback = e.feedbackValues(tIdx)'>0;
n.correctResponse = uint8(correctResponse(tIdx)>0)+1;

response = double(n.correctResponse);
response(n.feedback==0) = uint8(-1*(response(n.feedback==0)-2)+1);

galvoPosition = p.galvoCoords(abs(e.galvoPosValues(tIdx))',:);
galvoPosition(:,1) = galvoPosition(:,1).*sign(e.galvoPosValues(tIdx))';

n.response = uint8(response);
n.repeatNum = uint8(repeatNum');
n.rewardAvailable = e.rewardAvailableValues(tIdx)'>0;
n.reactionTime = single(reactionTime(tIdx));
n.clickDuration = single([v(tIdx).clickDuration]');
n.clickRate = single([v(tIdx).clickRate]');
n.audAmplitude = audAmplitude(tIdx,:);
n.visContrast = visContrast(tIdx,:);
n.audDiff = audDiff(tIdx,:);
n.visDiff = visDiff(tIdx,:);
n.audInitialAzimuth = audInitialAzimuth(tIdx);
n.visInitialAzimuth = visInitialAzimuth(tIdx);
n.visAltitude = [v(tIdx).visAltitude]';
n.visSigma = [v(tIdx).visSigma]';
n.galvoType = int8(e.galvoTypeValues(tIdx)');
n.galvoPosition = galvoPosition;
n.laserType = int8(e.laserTypeValues(tIdx)');
n.laserPower = single(e.laserPowerValues(tIdx)');
n.laserSession = sum(unique(n.laserType))+(n.laserPower*0);
n.laserOnOff = [e.galvoTTLTimes(tIdx)' e.galvoAndLaserEndTimes(tIdx)']-n.trialStart;


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
n.conditions = rightConditions + -1*leftConditions;

uniqueConditions = uniqueConditions(~ismember(uniqueConditions, [rightInitialConditions; leftInitialConditions], 'rows'),:);
uniqueConditions = [uniqueConditions; rightInitialConditions; leftInitialConditions];
uniqueConditionsIdx = [0*find(n.conditions==0,1); (1:size(rightInitialConditions,1))';  -1*(1:size(leftInitialConditions,1))'];

audTrial = (n.visContrast==0 | n.visInitialAzimuth==0) & (n.audAmplitude>0 & n.audInitialAzimuth~=0);
visTrial = (n.audAmplitude==0 | n.audInitialAzimuth==0) & (n.visContrast>0 & n.visInitialAzimuth~=0);
coherentTrial = sign(n.visInitialAzimuth.*n.audInitialAzimuth)>0 & n.audAmplitude>0 & n.visContrast>0;
coflictTrial = sign(n.visInitialAzimuth.*n.audInitialAzimuth)<0 & n.audAmplitude>0 & n.visContrast>0;
blankTrial = (n.audAmplitude == 0 | n.audInitialAzimuth == 0) & (n.visContrast == 0 | n.visInitialAzimuth == 0);
n.trialType = ~blankTrial.*(audTrial+visTrial*2+coherentTrial*3+coflictTrial*4);

r.rawWheelTimeValue = n.rawWheelTimeValue;
n = rmfield(n, 'rawWheelTimeValue');
n.uniqueConditions = uniqueConditions;
n.uniqueConditionsIdx = uniqueConditionsIdx;

p.maxRepeatIncorrect = max(p.maxRepeatIncorrect);
p.numberConditions = length(unique([audAmplitude, visContrast audInitialAzimuth visInitialAzimuth], 'rows'));
p.audPerformance = round(mean(n.feedback(n.trialType==1 & n.laserType == 0))*100);
p.visPerformance = round(mean(n.feedback(n.trialType==2 & n.laserType == 0))*100);
p.mulPerformance = round(mean(n.feedback(n.trialType==3 & n.laserType == 0))*100);
p.validTrials = sum(tIdx);

x.validTrials = tIdx;
newParams = p;
newBlock = n;
newRaw = r;

%% Check that all expected fields exist
blockFields = {'subject';'expDate';'sessionNum';'rigName';'rigType';'trialStart';'trialEnd';'stimPeriodStart';'closedLoopStart';'feedback'; ...
    'correctResponse';'response';'repeatNum';'reactionTime';'clickDuration';'clickRate';'audAmplitude';'visContrast';'audDiff';'visDiff';...
    'audInitialAzimuth';'visInitialAzimuth';'visAltitude';'visSigma';'galvoType';'galvoPosition';'laserType';'laserPower';'laserSession';...
    'laserOnOff';'conditions';'trialType';'uniqueConditions'; 'rewardAvailable'};
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
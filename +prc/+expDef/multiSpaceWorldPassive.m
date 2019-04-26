function x = multiSpaceWorldPassive(x)
%% A helper function for multisensoySpaceWorld experimental definition that produces standardised files with useful structures for further analysis.

%% Convert to shorter names for ease of use later
v = x.standardizedBlock.paramsValues;  %Parameter values at start of trial
e = x.standardizedBlock.events;        %Event structure
n = x.newBlock;                        %x.newBlock, already populated with subject, expDate, expNum, rigName, and expType
p = x.standardizedParams;              %Parameter values at start of entire session (includes multiple values for different conditions

p.visContrast = p.visContrast(1,:);
trialTimes = [e.newTrialTimes' e.endTrialTimes'];
p.audInitialAzimuth(p.audAmplitude == 0) = inf;       %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
p.visInitialAzimuth(p.visContrast == 0) = inf;        %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)

stimPeriodStart = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues == 1)';
%Process the raw data to be stored separately (because it is large). These are the times at which the visual and auditory stimuli turned on/off and
%all the wheel values. All are indexed by trial using the "indexByTrial" function, and times are relative to stimulus onset.
r.visStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues'], stimPeriodStart, [1 0]);
r.audStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues'], stimPeriodStart, [1 0]);
r.wheelTimeValue = prc.indexByTrial(trialTimes, n.rawWheelTimeValue(:,1), [n.rawWheelTimeValue(:,1) n.rawWheelTimeValue(:,2)], stimPeriodStart, [1 0]);
r.wheelTimeValue(cellfun(@isempty, r.wheelTimeValue)) = deal({[0 0]});
r.wheelTimeValue = cellfun(@(x) [x(:,1) x(:,2)-x(find([x(1:end-1,1);1]>0,1),2)], r.wheelTimeValue, 'uni', 0);
n = rmfield(n, 'rawWheelTimeValue');

n.trialStartEnd = [e.newTrialTimes' e.endTrialTimes'];
n.audAmplitude = [v.audAmplitude]';
n.audInitialAzimuth = [v.audInitialAzimuth]';     
n.audInitialAzimuth(n.audAmplitude == 0) = inf; 

n.visContrast = max([v.visContrast])';
n.visInitialAzimuth = [v.visInitialAzimuth]';
n.visInitialAzimuth(n.visContrast==0) = inf;
n.visAltitude = [v.visAltitude]';
n.visSigma = [v.visSigma]';
n.rewardTriggered = [v.feedback]';

%allConditions is all the conditions the mouse actually performed, where conditions can be completely defined by audAmplitude, visContrast, and the
%intial azimuth of aud and vis stimuli.
%Unique conditions is based on the parameter set, and includes all possible conditions that the mouse could have experienced. We repeat audAmplitude
%because in some cases it will only have one value, and in some cases it will have more;
allConditions = [n.audAmplitude n.visContrast n.audInitialAzimuth n.visInitialAzimuth];
uniqueConditions = unique([p.audAmplitude' p.visContrast' p.audInitialAzimuth' p.visInitialAzimuth'], 'rows');
%Create a set of unique conditions, where each row is a condition in the order: [zero conditions; right conditions; left conditions].
leftInitialConditions = uniqueConditions(uniqueConditions(:,end-1)< 0 | ((isinf(uniqueConditions(:,end-1)) | ~(uniqueConditions(:,end-1))) & uniqueConditions(:,end)<0),:);
rightInitialConditions = [leftInitialConditions(:,1:2) leftInitialConditions(:,end-1:end)*-1];
rightInitialConditions(isinf(rightInitialConditions)) = inf;
zeroConditions = uniqueConditions(~ismember(uniqueConditions, [leftInitialConditions; rightInitialConditions], 'rows'),:);
uniqueConditions = [zeroConditions; rightInitialConditions; leftInitialConditions];

%For each trial, find which row of leftInitialConditions or rightInitialConditions it belongs to, and check that no rows belong to both. We generate a
%condition index where each trial is positive if right, -ve if left, and the opposite of any conditionIdx is the inverse sign of that conditionIdx. We
%also create uniqueConditionReference which has the corresponding condition for each row of uniqueConditions. Finally, we create conditionRowIdx which (for
%every trial) inditcates which row of the uniqueConditions table that trial corresponds to.
[~, rightConditionsIdx] = ismember(allConditions, rightInitialConditions, 'rows');
[~, leftConditionsIdx] = ismember(allConditions, leftInitialConditions, 'rows');
if any(all([rightConditionsIdx~=0, leftConditionsIdx~=0],2)); error('Detect same condition as being Left and Right'); end
conditionLabel = rightConditionsIdx + -1*leftConditionsIdx;
if size(zeroConditions,1)>1
    [~, zeroConditionsIdx] = ismember(allConditions, zeroConditions, 'rows');
    zeroConditionsIdx(zeroConditionsIdx~=0) = zeroConditionsIdx(zeroConditionsIdx~=0)+max(rightConditionsIdx);
    conditionLabel = conditionLabel+zeroConditionsIdx;
end
uniqueConditionLabels = [unique(zeroConditionsIdx(zeroConditionsIdx~=0)); (1:size(rightInitialConditions,1))'; -1*(1:size(leftInitialConditions,1))'];
[~, conditionRowIdx] = ismember(conditionLabel, uniqueConditionLabels);

n.uniqueConditions = uniqueConditions;
n.uniqueConditionLabels = uniqueConditionLabels;
n.conditionLabelRow = [conditionLabel conditionRowIdx]; 

x.newParams = prc.chkThenRemoveFields(p, {'postQuiescentDelay';'laserOnsetDelays';'waveformType';'maxRepeatIncorrect';'galvoType';'laserPower';...#
    'laserTypeProportions';'laserDuration';'rewardTotal';});
x.newBlock = n;
x.newRaw = r;

%% Check that all expected fields exist
expectedBlockFields = {'subject';'expDate' ;'expNum';'expDef';'rigName';'expType';'trialStartEnd' ;'audAmplitude';'audInitialAzimuth' ;'visContrast';'visInitialAzimuth';...
    'visAltitude' ;'visSigma';'rewardTriggered' ;'uniqueConditions';'uniqueConditionLabels';'conditionLabelRow'};
expectedParamFields =  {'audAmplitude';'audInitialAzimuth';'backgroundNoiseAmplitude';'clickDuration';'clickRate';'closedLoopOnsetToneAmplitude';'feedback';'responseWindow';...
    'rewardSize';'visAltitude';'visContrast';'visInitialAzimuth';'visSigma';'numRepeats';'subject';'expDate';'expNum';'expDef';'rigName';...
    'expType';'totalTrials';'minutesOnRig';'laserSession'};

if any(~contains(fields(x.newBlock), expectedBlockFields))
    excessfields = fields(x.newBlock);
    excessfields = excessfields(~contains(excessfields, expectedBlockFields));
    fprintf('Field mistmatch in block file %s - \n', excessfields{:});
    keyboard;
elseif any(~contains(expectedBlockFields, fields(x.newBlock)))
    excessfields = expectedBlockFields(~contains(expectedBlockFields, fields(x.newBlock)));
    fprintf('Field mistmatch in block file %s - \n', excessfields{:});
    keyboard;
end

if any(~contains(fields(x.newParams), expectedParamFields))
    excessfields = fields(x.newParams);
    excessfields = excessfields(~contains(excessfields, expectedParamFields));
    fprintf('Field mistmatch in Param file %s - \n', excessfields{:});
    keyboard;
elseif any(~contains(expectedParamFields, fields(x.newParams)))
    excessfields = expectedParamFields(~contains(expectedParamFields, fields(x.newParams)));
    fprintf('Field mistmatch in Param file %s - \n', excessfields{:});
    keyboard;
end
end
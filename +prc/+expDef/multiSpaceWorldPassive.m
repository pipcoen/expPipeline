function x = multiSpaceWorldPassive(x)
%% A helper function for multisensoySpaceWorld experimental definition that produces standardised files with useful structures for further analysis.
%% Convert to shorter names for ease of use later
v = x.standardizedBlock.paramsValues;  %Parameter values at start of trial
e = x.standardizedBlock.events;        %Event structure
n = x.newBlock;                        %x.newBlock, already populated with subject, expDate, expNum, rigName, and expType
p = x.standardizedParams;              %Parameter values at start of entire session (includes multiple values for different conditions
eIdx = 1:length(e.endTrialTimes);

p.visContrast = p.visContrast(1,:);
trialTimes = [e.newTrialTimes(eIdx)' e.endTrialTimes(eIdx)'];
stimPeriodStart = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues == 1)';
p.audInitialAzimuth(p.audAmplitude == 0) = inf;       %Change case when audAmplitude was 0 to have infinite azimuth (an indication of no azimuth value)
p.visInitialAzimuth(p.visContrast == 0) = inf;        %Change case when visContrast was 0 to have infinite azimuth (an indication of no azimuth value)

%Process the raw data to be stored separately (because it is large). These are the times at which the visual and auditory stimuli turned on/off and
%all the wheel values. All are indexed by trial using the "indexByTrial" function, and times are relative to stimulus onset.
times2Subtract = [stimPeriodStart stimPeriodStart*0];
r.visStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues'], times2Subtract);
r.audStimOnOffTTimeValue = prc.indexByTrial(trialTimes, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues'], times2Subtract);
r.wheelTimeValue = prc.indexByTrial(trialTimes, n.rawWheelTimeValue(:,1), [n.rawWheelTimeValue(:,1) n.rawWheelTimeValue(:,2)], times2Subtract);
r.wheelTimeValue(cellfun(@isempty, r.wheelTimeValue)) = deal({[0 0]});
r.wheelTimeValue = cellfun(@(x) [x(:,1) x(:,2)-x(find([x(1:end-1,1);1]>0,1),2)], r.wheelTimeValue, 'uni', 0);
n = rmfield(n, 'rawWheelTimeValue');

audAmplitude = [v.audAmplitude]';
visContrast = max([v.visContrast])';
audInitialAzimuth = [v.audInitialAzimuth]';     
audInitialAzimuth(audAmplitude == 0) = inf; 
visInitialAzimuth = [v.visInitialAzimuth]';     
visInitialAzimuth(visContrast == 0) = inf; 
rewardClick = [v.feedback]'>0; 
closedLoopOnsetTone = [v.closedLoopOnsetToneAmplitude]'>0;

%allConditions is all the conditions the mouse actually performed, where conditions can be completely defined by audAmplitude, visContrast, and the
%intial azimuth of aud and vis stimuli.
%Unique conditions is based on the parameter set, and includes all possible conditions that the mouse could have experienced. We repeat audAmplitude
%because in some cases it will only have one value, and in some cases it will have more;
allConditions = [audAmplitude visContrast audInitialAzimuth visInitialAzimuth];
uniqueConditions = unique(allConditions(~(rewardClick | closedLoopOnsetTone),:), 'rows');
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
conditionLabel((rewardClick | closedLoopOnsetTone)) = nan;
uniqueConditionLabels = [0*find(zeroConditions,1); (1:size(rightInitialConditions,1))'; -1*(1:size(leftInitialConditions,1))'];
[~, conditionRowIdx] = ismember(conditionLabel, uniqueConditionLabels);
uniqueDiff = [uniqueConditions(:,3) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];

%Create a "trialType" field which is 0,1,2,3,4 for blank, auditory, visual, coherent, and incoherent trials.
trialClass.blank = (visContrast==0 | visInitialAzimuth==0) & (audAmplitude==0 | audInitialAzimuth==0);
trialClass.auditory = (visContrast==0 | visInitialAzimuth==0) & (audAmplitude>0 & audInitialAzimuth~=0);
trialClass.visual = (audAmplitude==0 | audInitialAzimuth==0) & (visContrast>0 & visInitialAzimuth~=0);
trialClass.coherent = sign(visInitialAzimuth.*audInitialAzimuth)>0 & audAmplitude>0 & visContrast>0;
trialClass.conflict = sign(visInitialAzimuth.*audInitialAzimuth)<0 & audAmplitude>0 & visContrast>0;
trialClass.rewardClick = rewardClick; 
trialClass.closedLoopOnsetTone = closedLoopOnsetTone;

[audDiff, visDiff] = deal(nan*conditionRowIdx); 
audDiff(conditionRowIdx~=0) = uniqueDiff(conditionRowIdx(conditionRowIdx~=0), 1);
visDiff(conditionRowIdx~=0) = uniqueDiff(conditionRowIdx(conditionRowIdx~=0), 2);

[grids.visValues, grids.audValues] = meshgrid(unique(uniqueDiff(:,2)), unique(uniqueDiff(:,1)));
[~, gridIdx] = ismember(uniqueDiff, [grids.audValues(:) grids.visValues(:)], 'rows');
grids.conditionLabels = nan*grids.visValues;
% grids.conditionLabels(gridIdx) = uniqueConditionLabels;

n.conditionParametersAV = uniqueDiff;
n.conditionLabels = uniqueConditionLabels;
n.timings.trialStartEnd = trialTimes;
n.trialClass = trialClass; 
n.timings.stimPeriodStart = stimPeriodStart;
n.stim.audAmplitude = audAmplitude;
n.stim.audInitialAzimuth = audInitialAzimuth;
n.stim.audDiff = audDiff;
n.stim.visContrast = visContrast;
n.stim.visInitialAzimuth = visInitialAzimuth;
n.stim.visDiff = visDiff;
n.stim.conditionLabel = conditionLabel; 
n.grids = grids;
n.params = x.oldParams;

x.newBlock = n;
x.newRaw = r;
end
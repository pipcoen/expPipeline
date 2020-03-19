function x = multiSpaceWorld(x)
%% A helper function for multisensoySpaceWorld experimental definition that produces standardised files with useful structures for further analysis.

% INPUTS
% x---------------------Structure form the "convertExpFiles" function that contains all the file information

% OUTPUTS
% "x.newBlock" is the new, compact, and restructured block file with the following fields:
    %.subject---------------Name of the mouse
    %.expDate---------------Date that the experiment was recorded
    %.expNum----------------Session number for experiment
    %.rigName---------------Name of the rig where the experiment took place
    %.expType---------------Type experiment ('training, inactivation, or ephys')
    %.expDef----------------Name of the expDef file
    %.performanceAVM--------Performance of mouse (%) on auditory, visual, and multisensory trials
    %.conditionParametersAV-[visDiff, audDiff], one row for each condition in the parameter set for the session
    %.conditionLabels-------Integer labels, one for each condition
    
    %.trialClass------------Structure with logical fields to identifying each class of trial
        %.blank----------------Trials with auditory in the center and visual contrast is zero
        %.auditory-------------Trials with auditory on left/right and visual contrast is zero
        %.visual---------------Trials with auditory in the center and visual contrast is not zero
        %.coherent-------------Trials where the auditory and visual stimuli agree (and are non-zero)
        %.conflict-------------Trials with the auditory and visual stimuli disagree (and are non-zero)
        
    %.timings---------------Structure containing timings for events on each trial
        %.trialStartEnd--------[start end] times for whole trial
        %.stimPeriodStart------Time of stimulus onset (aud and vis start at the same time)
        %.closedLoopStart------Time of closed loop initiation (typically 500ms after stimulus onset)
   
    %.stim------------------Structure containing information about the stimulus presented on each trial
        %.audAmplitude---------Aud amplitude (arbitrary number really, only ~=0 matters)
        %.audInitialAzimuth----Initial azimuthal location of auditory stimulus (+/- is right/left of mouse). Inf if not present.
        %.audDiff--------------Difference in initial azimuthal location (identical to audInitialAzimuth in most cases)
        %.visContrast----------Absolute visual contrast (+ for left and right). Ranges from 0 to 1
        %.visInitialAzimuth----Initial azimuthal location of visual stimulus (+/- is right/left of mouse). Inf if not present.
        %.visDiff--------------Difference in left/right contrast (+/- is right/left of mouse)
        %.conditionLabel-------The integer label fpr the condition being presented.
        
    %.inactivation----------Structure containing information about the inactivation paramters for each trial
        %.laserType------------Laser type, which can be off (0), unilateral (1), or bilateral (2)
        %.laserPower-----------Laser power (represents mW in more recent experiments)
        %.galvoType------------Galvo movement type, which can be stationary (1) or moving between two sites (2) while the laser is on
        %.laserOnsetDelay------Delay between laser onset and stimulus onset
        %.galvoPosition--------Galvo coordinates used on each trial ([LM axis, AP axis], +/- is right/left hemisphere)
        %.laserOnOff-----------[on off] times for the laser, relative to trial start.

    %.outcome---------------Structure containing information about the outcome of each trial
        %.timeToWheelMove------Time between the stimulus onset and significant wheel movement
        %.responseMade---------[on off] times for the laser, relative to trial start.       
        
    %.params----------------The parameters from signals. Basically never used, but saved for completion.
        

% "x.newRaw" is a structure comprising potentially useful raw data (such as wheel movement and timeline data) which is not used for a lot of analyses and
% so should only be loaded if necessary (as it is large). Fields are:
    %.visStimOnOffTTimeValue---------cell array, each contatining an [nx2] vector of [time(s), value(on/off is 1/0)] for the visual stimulus
    %.audStimOnOffTTimeValue---------cell array, each contatining an [nx2] vector of [time(s), value(on/off is 1/0)] for the auditory stimulus
    %.wheelTimeValue-----------------wheel [times(s), position(deg)] over the course of the entire session
    %.visAzimuthTimeValue------------cell array, each contatining an [nx2] vector of [time(s), value(deg)] for the visual stimulus azimuth
    %.audAzimuthTimeValue------------cell array, each contatining an [nx2] vector of [time(s), value(deg)] for the auditory stimulus azimuth

%% Convert to shorter names for ease of use later
v = x.standardizedBlock.paramsValues;  %Parameter values at start of trial
e = x.standardizedBlock.events;        %Event structure
n = x.newBlock;                        %x.newBlock, already populated with subject, expDate, expNum, rigName, and expType
p = x.standardizedParams;              %Parameter values at start of entire session (includes multiple values for different conditions
vIdx = e.repeatNumValues(1:length(x.standardizedBlock.events.endTrialTimes))==1;            %Indices of valid trials (0 for repeats)

%% Trim trials if there are more than 100 total trials (if not, the mouse was likely still learning)
if sum(vIdx) > 150
    %We invalidate the first 5 and last 5 correct trials for each session.
    vIdx = double(vIdx);
    stimStartTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1);
    quickResponses = (e.feedbackTimes(1:length(vIdx)) - stimStartTimes(1:length(vIdx)))<1.5;
    vIdx(1:max(find(vIdx==1 & e.responseTypeValues(1:length(vIdx))~=0 & quickResponses, 5, 'first'))) = -1;
    vIdx(min(find(vIdx==1 & e.responseTypeValues(1:length(vIdx))~=0 & quickResponses, 5, 'last')):end) = -1;
    
    %Invalidate trials where laser "trasitionTimes" are more than 90% of the time until the TTL pulse that activates the laser. Otherwise, cannot
    %be confident that the laser was ready to receive the pulse.
    if strcmp(n.expType, 'inactivation')
        trasitionTimes = e.laserInitialisationTimes(1:length(vIdx))-e.newTrialTimes(1:length(vIdx));
        timeToTTL = (e.galvoTTLTimes(1:length(vIdx))-e.newTrialTimes(1:length(vIdx)))*0.9;
        vIdx(trasitionTimes(:)>timeToTTL(:) & vIdx(:)==1)=-1;
    end
    vIdx = vIdx>0;
end

%% The number of repeats and timeouts for each trial type presented
%Invalidate trials that are repeats following an incorrect choice (because the mouse knows which way to go based on the incorrect choice) and trials
%where there were multiple repeats because of a timeout (i.e. only the initial timeout response is counted as valid)
for i = find(vIdx)
    if i < length(vIdx) && e.responseTypeValues(i) == 0
        nextResponse = min([i+find(e.responseTypeValues(i+1:length(vIdx))~=0 | e.repeatNumValues(i+1:length(vIdx))==1,1), length(vIdx)+1]);
        if e.repeatNumValues(nextResponse)==1 || nextResponse >= length(vIdx); continue; end
        vIdx(nextResponse) = 1;
    end
end

%% Extract meaningful data from the block file
%eIdx is just an logical for all trials that ended (if the experiment ends mid-trial, there may be an extra index for some events)
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

%Get trial start/end times, stim start times, closed loop start times, feedback times, etc.
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

%Calculate an approximate time to the first wheel movement. This is different from the "timeToFeedback" in that it is based on wheel movement, rather
%than the time when the threshold was reached. WheelMove is an interpolation of the wheel movement (to get it's nearest position at every ms).
%wheelToThresh uses the difference between the wheel position at the closed loop start and threshold to calculate the change in wheel value that
%represents a response (this can be very different for different rotary encoders). timeToWheelMove is then the time at which wheelMove exceeds 25% of
%wheelToThresh
%%
time2closedLoopResponded = closedLoopStart(~timeOuts)-stimPeriodStart(~timeOuts);
timeToFeedbackResponded = timeToFeedback(~timeOuts);
wheelPos = cellfun(@(x) [(-1:0.001:x(end,1))', (interp1(x(diff(x(:,1))~=0,1), x(diff(x(:,1))~=0,2), -1:0.001:x(end,1), 'nearest', 'extrap'))'], r.wheelTimeValue(~timeOuts), 'uni', 0);

wheelToThresh = arrayfun(@(x,y,z) x{1}(x{1}(:,1)>y & x{1}(:,1)<z,2), wheelPos, time2closedLoopResponded, timeToFeedbackResponded, 'uni', 0);
wheelToThresh = nanmedian(abs(cellfun(@(x) x(end)-x(1), wheelToThresh(~cellfun(@isempty, wheelToThresh)))));
time2Thresh = (cellfun(@(x) x(find(abs(x(:,2))>wheelToThresh/4 & x(:,1)>0,1),1), wheelPos, 'uni', 0));

smoothWindow = 31;
velThresh  = max([5 wheelToThresh/25]);

wheelSmooth = cellfun(@(x) [x(:,1) smooth(x(:,2), smoothWindow)], wheelPos, 'uni', 0);
wheelVel = cellfun(@(x) [x(:,1) [0; diff(x(:,2))]*1000], wheelSmooth, 'uni', 0);
velAboveMinThresh = cellfun(@(x) [x(:,1) convn(abs(x(:,2)')>velThresh, ones(1,smoothWindow), 'same')'], wheelVel, 'uni', 0);
velOnsetTimes = cellfun(@(x) x(strfind(x(:,2)'>=smoothWindow, [0,1]),1), velAboveMinThresh, 'uni', 0);
velOnsetTimes = (cellfun(@(x,y) max([nan x(find(x>0 & x<y, 1, 'last'))]), velOnsetTimes, time2Thresh, 'uni', 0));
timeToWheelMove = feedbackValues;
timeToWheelMove(~timeOuts) = cell2mat(velOnsetTimes);

%Get the response the mouse made on each trial based on the correct response and then taking the opposite for incorrect trials. NOTE: this will not
%work for a task with more than two response options.
responseMade = double(correctResponse).*~timeOuts;
responseMade(feedbackValues<0) = -1*(responseMade(feedbackValues<0));

%allConditions is all the conditions the mouse actually performed, where conditions can be completely defined by audAmplitude, visContrast, and the
%initial azimuth of aud and vis stimuli.
%Unique conditions is based on the parameter set, and includes all possible conditions that the mouse could have experienced. We repeat audAmplitude
%because in some cases it will only have one value, and in some cases it will have more;
allConditions = [audAmplitude visContrast audInitialAzimuth visInitialAzimuth];
uniqueConditions = unique([p.audAmplitude' p.visContrast' p.audInitialAzimuth' p.visInitialAzimuth'], 'rows');

%Create a set of unique conditions, where each row is a condition in the order: [zero conditions; right visual conditions; left visual conditions].
leftInitialConditions = uniqueConditions(uniqueConditions(:,end-1)< 0 | ((isinf(uniqueConditions(:,end-1)) | ~(uniqueConditions(:,end-1))) & uniqueConditions(:,end)<0),:);
if size(leftInitialConditions,1)~=floor(size(uniqueConditions,1)/2); warning('Why are half conditions not left conditions?'); end
zeroConditions = uniqueConditions(all([any(uniqueConditions(:,[1,3])==0,2) any(uniqueConditions(:,[2,4])==0,2)],2),:);
rightInitialConditions = [leftInitialConditions(:,1:2) leftInitialConditions(:,end-1:end)*-1];
rightInitialConditions(isinf(rightInitialConditions)) = inf;
rightInitialConditions = [rightInitialConditions; uniqueConditions(~ismember(uniqueConditions, [zeroConditions; leftInitialConditions; rightInitialConditions], 'rows'),:)];
uniqueConditions = [zeroConditions; rightInitialConditions; leftInitialConditions];

%For each trial, find which row of leftInitialConditions or rightInitialConditions it belongs to, and check that no rows belong to both. We generate a
%condition index where each trial is positive if right, -ve if left, and the opposite of any conditionIdx is the inverse sign of that conditionIdx. We
%also create uniqueConditionReference which has the corresponding condition for each row of uniqueConditions. Finally, we create conditionRowIdx which
%(for every trial) inditcates which row of the uniqueConditions table that trial corresponds to.
[~, rightConditionsIdx] = ismember(allConditions, rightInitialConditions, 'rows');
[~, leftConditionsIdx] = ismember(allConditions, leftInitialConditions, 'rows');
if any(all([rightConditionsIdx~=0, leftConditionsIdx~=0],2)); error('Detect same condition as being Left and Right'); end
conditionLabel = rightConditionsIdx + -1*leftConditionsIdx;
uniqueConditionLabels = [0*find(zeroConditions,1); (1:size(rightInitialConditions,1))'; -1*(1:size(leftInitialConditions,1))'];
uniqueDiff = [uniqueConditions(:,3) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];
[~, conditionRowIdx] = ismember(conditionLabel, uniqueConditionLabels);

%Create a "logical" for each trial type (blank, auditory, visual, coherent, and incoherent trials)
trialClass.blank = (visContrast==0 | visInitialAzimuth==0) & (audAmplitude==0 | audInitialAzimuth==0);
trialClass.auditory = (visContrast==0 | visInitialAzimuth==0) & (audAmplitude>0 & audInitialAzimuth~=0);
trialClass.visual = (audAmplitude==0 | audInitialAzimuth==0) & (visContrast>0 & visInitialAzimuth~=0);
trialClass.coherent = sign(visInitialAzimuth.*audInitialAzimuth)>0 & audAmplitude>0 & visContrast>0;
trialClass.conflict = sign(visInitialAzimuth.*audInitialAzimuth)<0 & audAmplitude>0 & visContrast>0;

audDiff = uniqueDiff(conditionRowIdx, 1);
visDiff = uniqueDiff(conditionRowIdx, 2);

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
n.params = x.oldParams;

x.newBlock = n;
%Remove fields of the raw data where the trial lasts longer than 5s. The data from there trials aren't likely to be useful (since the mouse stayed
%still for a long time) and can take up a lot of space in the .mat file.
for i = fields(r)'; newRaw.(i{1}) = r.(i{1}); newRaw.(i{1})(x.newBlock.outcome.timeToFeedback>5) = {[]}; end
x.newRaw = r;
end
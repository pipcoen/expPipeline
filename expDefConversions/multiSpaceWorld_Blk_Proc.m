function [newBlock, newParams, newRaw] = multiSpaceWorld_Blk_Proc(x)
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

p.galvoCoords = e.galvoCoordsValues(:,1:2);
correctResponse = [v.correctResponse]';  %Correct Response
audInitialAzimuth = [v.audInitialAzimuth]'; audInitialAzimuth(audAmplitude==0) = inf;
visInitialAzimuth = [v.visInitialAzimuth]'; visInitialAzimuth(visContrast==0) = inf;

audDiff = audAmplitude.*(audInitialAzimuth<0)-audAmplitude.*(audInitialAzimuth>0);
visDiff = visContrast.*(audInitialAzimuth<0)-visContrast.*(audInitialAzimuth>0);

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
n.laserOnOff = [e.galvoTTLTimes(tIdx)' e.galvoAndLaserEndTimes(tIdx)']-n.trialStart;


allConditions = [n.audAmplitude n.visContrast n.audInitialAzimuth n.visInitialAzimuth];
uniqueConditions = unique(allConditions, 'rows');
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

p.visContrast = unique(n.visContrast)';
p.audAmplitude = unique(n.audAmplitude)';
p.maxRepeatIncorrect = max(p.maxRepeatIncorrect);
p.visInitialAzimuth = unique(abs(n.visInitialAzimuth))';
p.audInitialAzimuth = unique(abs(n.audInitialAzimuth))';
p.numberConditions = length(unique([audAmplitude, visContrast audInitialAzimuth visInitialAzimuth], 'rows'));
p.audPerformance = round(mean(n.feedback(n.trialType==1))*100);
p.visPerformance = round(mean(n.feedback(n.trialType==2))*100);
p.mulPerformance = round(mean(n.feedback(n.trialType==3))*100);
p.validTrials = sum(tIdx);

x.validTrials = tIdx;
newParams = p;
newBlock = n;
newRaw = r;
end
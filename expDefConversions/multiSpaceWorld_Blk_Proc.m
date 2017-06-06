function [newBlock, newParams] = multiSpaceWorld_Blk_Proc(x)
v = x.standardizedBlock.paramsValues;
e = x.standardizedBlock.events;
n = x.newBlock;
p = x.standardizedParams;

correctResponse = [v.correctResponse]';  %Correct Response
audInitialAzimuth = [v.audInitialAzimuth]';
visInitialAzimuth = [v.visInitialAzimuth]';

audAmplitude = single(cell2mat({v.audAmplitude}'));
audLeftRight = [audAmplitude.*(audInitialAzimuth<0) audAmplitude.*(audInitialAzimuth>0)];

visContrast = single(cell2mat({v.visContrast}')); 
visLeftRight = [visContrast.*(visInitialAzimuth<0) visContrast.*(visInitialAzimuth>0)];

n.visStimOnOff = indexByTrial(n, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues']);
n.audStimOnOff = indexByTrial(n, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues']);

n.wheelTimeValue = indexByTrial(n, n.rawWheelTimeValue(:,1), n.rawWheelTimeValue);
n.visAzimuthTimeValue = indexByTrial(n, e.visAzimuthTimes', [e.visAzimuthTimes' e.visAzimuthValues']);
n.audAzimuthTimeValue = indexByTrial(n, e.audAzimuthTimes', [e.audAzimuthTimes' e.audAzimuthValues']);

closedLoopStart = e.closedLoopOnOffTimes(e.closedLoopOnOffValues == 1)';
reactionTime = e.feedbackTimes(x.validTrials)'-closedLoopStart(x.validTrials);

closedLoopStart = indexByTrial(n, closedLoopStart, closedLoopStart, 1);
n.closedLoopStart = single(cell2mat(closedLoopStart));

n.feedback = e.feedbackValues(x.validTrials)'>0; 
n.correctResponse = uint8(correctResponse(x.validTrials)>0)+1;

response = n.correctResponse;
response(n.feedback==0) = -1*(response(n.feedback==0)-2)+1;

n.response = uint8(response);
n.repeatNum = uint8(x.repeatNum');
n.reactionTime = single(reactionTime);
n.clickDuration = single([v(x.validTrials).clickDuration]');
n.clickRate = single([v(x.validTrials).clickRate]');
n.audAmplitude = audAmplitude(x.validTrials,:);
n.visContrast = visContrast(x.validTrials,:);
n.audLeftRight = audLeftRight(x.validTrials,:);
n.visLeftRight = visLeftRight(x.validTrials,:);
n.audInitialAzimuth = [v(x.validTrials).audInitialAzimuth]';
n.visInitialAzimuth = [v(x.validTrials).visInitialAzimuth]';
n.visAltitude = [v(x.validTrials).visAltitude]';
n.visSigma = [v(x.validTrials).visSigma]';

allConditions = [n.audLeftRight n.visLeftRight n.audInitialAzimuth.*(n.audAmplitude>0) n.visInitialAzimuth.*(n.visContrast>0)];
n.uniqueConditions = unique(allConditions, 'rows');
[~, n.conditions] = ismember(allConditions, n.uniqueConditions, 'rows');

audTrial = ~any(n.visContrast,2);
visTrial = ~any(n.audAmplitude,2);
coherentTrial = sign(n.visInitialAzimuth.*n.audInitialAzimuth)>0 & any(n.audAmplitude,2) & any(n.visContrast,2);
coflictTrial = sign(n.visInitialAzimuth.*n.audInitialAzimuth)<0 & any(n.audAmplitude,2) & any(n.visContrast,2);
centeredAudTrial = n.audAmplitude>0 & n.audInitialAzimuth==0;
blankTrial = ~any(n.audAmplitude,2)&~any(n.visContrast,2)*2;
n.trialType = ~blankTrial.*(audTrial+visTrial*2+coherentTrial*3+coflictTrial*4+centeredAudTrial*5);


p.visContrast = unique(n.visContrast)';
p.audAmplitude = unique(n.audAmplitude)';
p.maxRepeatIncorrect = max(p.maxRepeatIncorrect);
p.audInitialAzimuth = unique(abs(n.visInitialAzimuth))';
p.visInitialAzimuth = unique(abs(n.audInitialAzimuth))';
p.numberConditions = length(unique([audAmplitude, visContrast audInitialAzimuth visInitialAzimuth], 'rows'));
p.audPerformance = round(mean(n.feedback(n.trialType==1))*100);
p.visPerformance = round(mean(n.feedback(n.trialType==2))*100);
p.mulPerformance = round(mean(n.feedback(n.trialType==3))*100);

newParams = p;
newBlock = n;
end
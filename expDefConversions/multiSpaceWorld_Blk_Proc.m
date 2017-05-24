function [blk, prm] = multiSpaceWorld_Blk_Proc(x, b, blk, prm)
v = b.paramsValues;
e = b.events;

correctResponse = [v.correctResponse]';  %Correct Response
audInitialAzimuth = [v.audInitialAzimuth]';
visInitialAzimuth = [v.visInitialAzimuth]';

audAmplitude = single(cell2mat({v.audAmplitude}'));
audLeftRight = [audAmplitude.*(audInitialAzimuth<0) audAmplitude.*(audInitialAzimuth>0)];

visContrast = single(cell2mat({v.visContrast}')); 
visLeftRight = [visContrast.*(visInitialAzimuth<0) visContrast.*visInitialAzimuth>0];

blk.visStimOnOff = indexByTrial(blk, e.visStimOnOffTimes', [e.visStimOnOffTimes' e.visStimOnOffValues']);
blk.audStimOnOff = indexByTrial(blk, e.audStimOnOffTimes', [e.audStimOnOffTimes' e.audStimOnOffValues']);

blk.wheelTimeValue = indexByTrial(blk, blk.rawWheelTimeValue(:,1), blk.rawWheelTimeValue);
blk.visAzimuthTimeValue = indexByTrial(blk, e.visAzimuthTimes', [e.visAzimuthTimes' e.visAzimuthValues']);
blk.audAzimuthTimeValue = indexByTrial(blk, e.audAzimuthTimes', [e.audAzimuthTimes' e.audAzimuthValues']);

closedLoopStart = e.closedLoopOnOffTimes(e.closedLoopOnOffValues == 1)';
reactionTime = e.feedbackTimes(x.validTrials)'-closedLoopStart(x.validTrials);

closedLoopStart = indexByTrial(blk, closedLoopStart, closedLoopStart, 1);
blk.closedLoopStart = single(cell2mat(closedLoopStart));

blk.feedback = e.feedbackValues(x.validTrials)'>0; 
blk.correctResponse = uint8(correctResponse(x.validTrials)>0)+1;

response = blk.correctResponse;
response(blk.feedback==0) = -1*(response(blk.feedback==0)-2)+1;

blk.response = uint8(response);
blk.repeatNum = uint8(x.repeatNum');
blk.reactionTime = single(reactionTime);
blk.clickDuration = single([v(x.validTrials).clickDuration]');
blk.clickRate = single([v(x.validTrials).clickRate]');
blk.audAmplitude = audAmplitude(x.validTrials,:);
blk.visContrast = visContrast(x.validTrials,:);
blk.audLeftRight = audLeftRight(x.validTrials,:);
blk.visLeftRight = visLeftRight(x.validTrials,:);
blk.audInitialAzimuth = [v(x.validTrials).audInitialAzimuth]';
blk.visInitialAzimuth = [v(x.validTrials).visInitialAzimuth]';
blk.visAltitude = [v(x.validTrials).visAltitude]';
blk.visSigma = [v(x.validTrials).visSigma]';

audTrial = ~any(blk.visContrast,2);
visTrial = ~any(blk.audAmplitude,2);
coherentTrial = sign(blk.visInitialAzimuth.*blk.audInitialAzimuth)>0 & any(blk.audAmplitude,2) & any(blk.visContrast,2);
coflictTrial = sign(blk.visInitialAzimuth.*blk.audInitialAzimuth)<0 & any(blk.audAmplitude,2) & any(blk.visContrast,2);
centeredAudTrial = blk.audAmplitude>0 & blk.audInitialAzimuth==0;
blankTrial = ~any(blk.audAmplitude,2)&~any(blk.visContrast,2)*2;
blk.trialType = ~blankTrial.*(audTrial+visTrial*2+coherentTrial*3+coflictTrial*4+centeredAudTrial*5);


prm.visContrast = unique(blk.visContrast)';
prm.audAmplitude = unique(blk.audAmplitude)';
prm.maxRepeatIncorrect = max(prm.maxRepeatIncorrect);
prm.audInitialAzimuth = unique(abs(blk.visInitialAzimuth))';
prm.visInitialAzimuth = unique(abs(blk.audInitialAzimuth))';
prm.numberConditions = length(unique([audAmplitude, visContrast audInitialAzimuth visInitialAzimuth], 'rows'));
prm.audPerformance = round(mean(blk.feedback(blk.trialType==1))*100);
prm.visPerformance = round(mean(blk.feedback(blk.trialType==2))*100);
prm.mulPerformance = round(mean(blk.feedback(blk.trialType==3))*100);
end
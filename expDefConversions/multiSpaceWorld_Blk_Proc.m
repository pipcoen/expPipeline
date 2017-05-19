function [blk, prm] = multiSpaceWorld_Blk_Proc(x, b, blk, prm)
b.paramsValues = b.paramsValues;
correctResponse = [b.paramsValues.correctResponse]';  %Correct Response
audAmplitude = single(cell2mat({b.paramsValues.a
audInitialAzimuth = [b.paramsValues.audInitialAzimuth]';
visInitialAzimuth = [b.paramsValues.visInitialAzimuth]';


audLeftRight = [audAmplitude.*double(audInitialAzimuth<0) audAmplitude.*double(audInitialAzimuth>0)];
visContrast = single(cell2mat({b.paramsValues.visualContrast}')); 
visLeftRight = [visContrast.*double(visInitialAzimuth<0) visContrast.*double(visInitialAzimuth>0)];
vASg = [b.paramsValues.visualAltitudeSigma]';

stmV = sigOnOffTimes(b.events.stmVValues, b.events.stmVTimes);
blk.stmV = indexByTrial(blk, stmV(:,1), stmV, [1,1]);
stmA = sigOnOffTimes(b.events.stmAValues, b.events.stmATimes);
blk.stmA = indexByTrial(blk, stmA(:,1), stmA(:,1), 1);

whTV = indexByTrial(blk, blk.wrTV(:,1), blk.wrTV, [1,2]);
vsTV = [b.events.vPosTimes', b.events.vPosValues'];
vsTV = indexByTrial(blk, vsTV(:,1), vsTV, [1 0]);
asTV = [b.events.aPosTimes', b.events.aPosValues'];
asTV = indexByTrial(blk, asTV(:,1), asTV, [1 0]);
intO = b.events.intOTimes(x.validTrials)'-blk.sSrt;

fBck = b.events.fBckValues(x.validTrials)'>0; 
correctResponse = single(correctResponse(x.validTrials)>0)+1;
resp = correctResponse;
resp(fBck==0) = -1*(resp(fBck==0)-2)+1;

blk.correctResponse = correctResponse;
blk.resp = resp;
blk.rNum = uint8(x.rNum');
blk.fBck = fBck;
blk.rTim = single(b.events.fBckTimes(x.validTrials)'-b.events.intOTimes(x.validTrials)');
blk.clickDuration = single([b.paramsValues(x.validTrials).clickDuration]');
blk.clickRate = single([b.paramsValues(x.validTrials).clickRate]');
blk.audAmplitude = audAmplitude(x.validTrials,:);
blk.vCon = vCon(x.validTrials,:);
blk.aLeR = aLeR(x.validTrials,:);
blk.vLeR = vLeR(x.validTrials,:);
blk.audInitialAzimuth = [b.paramsValues(x.validTrials).audInitialAzimuth]';
blk.visInitialAzimuth = [b.paramsValues(x.validTrials).visInitialAzimuth]';
blk.vASg = single(vASg(x.validTrials,:));
blk.intO = single(intO);
blk.whTV = whTV;
blk.asTV = asTV;
blk.vsTV = vsTV;

audT = ~any(blk.vCon,2);
visT = ~any(blk.audAmplitude,2);
cohT = sign(blk.vIni.*blk.aIni)>0 & any(blk.audAmplitude,2) & any(blk.vCon,2);
conT = sign(blk.vIni.*blk.aIni)<0 & any(blk.audAmplitude,2) & any(blk.vCon,2);
cenT = blk.audAmplitude>0 & blk.aIni==0;
blnT = ~any(blk.audAmplitude,2)&~any(blk.vCon,2)*2;
blk.tTyp = ~blnT.*(audT+visT*2+cohT*3+conT*4+cenT*5);


prm.visualContrast = unique(prm.visualContrast)';
prm.clickDuration = prm.clickDurRate(1);
prm.clickRate = prm.clickDurRate(2);
prm.preStimQuiDuration = prm.preStimQuiRangeThr(1);
prm.audAmplitude = unique(prm.audAmplitude)';
prm.audioAzimuth = unique(abs(prm.audioAzimuth))';
prm.visualAzimuth = unique(abs(prm.visualAzimuth))';
prm.numberConditions = length(unique([audAmplitude, vCon], 'rows'));
prm.visualPerformance = round(mean(fBck(blk.tTyp==1))*100);
prm.audioPerformance = round(mean(fBck(blk.tTyp==2))*100);
prm.multiPerformance = round(mean(fBck(blk.tTyp==3))*100);

f2Re = {'interTrialDelay';'numRepeats';'interactPunishDelays'; ...
   'clickDurRate'; 'visualAltitudeSigma'; 'preStimQuiRangeThr'};
prm = chkThenRemoveFields(prm, f2Re);
end
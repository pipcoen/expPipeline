function [blk, prm] = multiSpaceWorld_Blk_Proc(x, b, blk, prm)
pVals = b.paramsValues;
cRes = [pVals.correctResponse]';  %Correct Response
ckDR = [pVals.clickDurRate]';     %Click Duration and Rate
aAmp = single(cell2mat({pVals.audioAmplitude}'));   %AudioAmplitude
aIni = [b.paramsValues.audioAzimuth]';
vIni = [b.paramsValues.visualAzimuth]';
aLeR = [aAmp.*double(aIni<0) aAmp.*double(aIni>0)];
vCon = single(cell2mat({pVals.visualContrast}')); 
vLeR = [vCon.*double(vIni<0) vCon.*double(vIni>0)];
vASg = [pVals.visualAltitudeSigma]';

stmV = sigOnOffTimes(b.events.stmVValues, b.events.stmVTimes);
blk.stmV = indexByTrial(blk, stmV(:,1), stmV, [1,1]);
stmA = sigOnOffTimes(b.events.stmAValues, b.events.stmATimes);
blk.stmA = indexByTrial(blk, stmA(:,1), stmA(:,1), 1);

whTV = indexByTrial(blk, blk.wrTV(:,1), blk.wrTV, [1,2]);
vsTV = [b.events.vPosTimes', b.events.vPosValues'];
vsTV = indexByTrial(blk, vsTV(:,1), vsTV, [1 0]);
asTV = [b.events.aPosTimes', b.events.aPosValues'];
asTV = indexByTrial(blk, asTV(:,1), asTV, [1 0]);
intO = b.events.intOTimes(x.vTri)'-blk.sSrt;

fBck = b.events.fBckValues(x.vTri)'>0; 
cRes = single(cRes(x.vTri)>0)+1;
resp = cRes;
resp(fBck==0) = -1*(resp(fBck==0)-2)+1;

blk.cRes = cRes;
blk.resp = resp;
blk.rNum = uint8(x.rNum');
blk.fBck = fBck;
blk.rTim = single(b.events.fBckTimes(x.vTri)'-b.events.intOTimes(x.vTri)');
blk.ckDR = single(ckDR(x.vTri,:));
blk.aAmp = aAmp(x.vTri,:);
blk.vCon = vCon(x.vTri,:);
blk.aLeR = aLeR(x.vTri,:);
blk.vLeR = vLeR(x.vTri,:);
blk.aIni = aIni(x.vTri);
blk.vIni = vIni(x.vTri);
blk.vASg = single(vASg(x.vTri,:));
blk.intO = single(intO);
blk.whTV = whTV;
blk.asTV = asTV;
blk.vsTV = vsTV;

audT = ~any(blk.vCon,2);
visT = ~any(blk.aAmp,2);
cohT = sign(blk.vIni.*blk.aIni)>0 & any(blk.aAmp,2) & any(blk.vCon,2);
conT = sign(blk.vIni.*blk.aIni)<0 & any(blk.aAmp,2) & any(blk.vCon,2);
cenT = blk.aAmp>0 & blk.aIni==0;
blnT = ~any(blk.aAmp,2)&~any(blk.vCon,2)*2;
blk.tTyp = ~blnT.*(audT+visT*2+cohT*3+conT*4+cenT*5);


prm.visualContrast = unique(prm.visualContrast)';
prm.clickDuration = prm.clickDurRate(1);
prm.clickRate = prm.clickDurRate(2);
prm.preStimQuiDuration = prm.preStimQuiRangeThr(1);
prm.audioAmplitude = unique(prm.audioAmplitude)';
prm.audioAzimuth = unique(abs(prm.audioAzimuth))';
prm.visualAzimuth = unique(abs(prm.visualAzimuth))';
prm.numberConditions = length(unique([aAmp, vCon], 'rows'));
prm.visualPerformance = round(mean(fBck(blk.tTyp==1))*100);
prm.audioPerformance = round(mean(fBck(blk.tTyp==2))*100);
prm.multiPerformance = round(mean(fBck(blk.tTyp==3))*100);

f2Re = {'interTrialDelay';'numRepeats';'interactPunishDelays'; ...
   'clickDurRate'; 'visualAltitudeSigma'; 'preStimQuiRangeThr'};
prm = chkThenRemoveFields(prm, f2Re);
end
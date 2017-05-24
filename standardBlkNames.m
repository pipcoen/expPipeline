function [b, p] = standardBlkNames(b, p)
% All trials are valid trials if there are no repeats.
nTri = length(b.paramsValues);
e = b.events;
v = b.paramsValues;

f2Re = {'audDevIdx';'audSampleRate';'numAudChannels';'type'; ...
    'services'; 'defFunction'; 'experimentIdx'; 'correctResponse';...
    'servicesDescription'; 'clickAmpDurRate'; 'vStimAltitude'; ...
    'stimulusAzimuth'; 'audVisAzimuth'; 'vStimSigma'; ...
    'interactPunishDelays'; 'stimulusDurRep'; 'noiseBurstAmpDur';...
    'rewardDurSize'; 'interactSigOnDurAmp'; 'audVisThreshold'; ...
    'stimulusContrast'; 'maxRetryIfIncorrect'; 'backNoiseAmp'; ...
    'preStimQuiRangeThr'; 'interTrialDelay';......
    ...
    ...
    'sPosTimes'; 'sPosValues'; 'stimStartTimes'; 'stimStartValues';...
    'sSrtTimes'; 'sSrtValues'; 'stmVValues'; 'stmVTimes'; ...
    'stmAValues'; 'stmATimes'; 'intOValues'; 'intOTimes'; ...
    'fBckValues'; 'fBckTimes'; 'rTotValues'; 'rTotTimes';...
    'visInitialAzimuthValues'; 'visInitialAzimuthTimes';...
    'audInitialAzimuthValues'; 'audInitialAzimuthTimes';...
    'preStimQuiescentDurationValues'; 'preStimQuiescentDurationTimes'};

if isfield(e, 'fBckTimes'); e.feedbackTimes = e.fBckTimes; e.feedbackValues = e.fBckValues; end
if isfield(e, 'stimStartTimes'); e.sSrtTimes = e.stimStartTimes; end
if isfield(e, 'stimStartTimes'); e.sSrtTimes = e.stimStartTimes; end
if isfield(p, 'backNoiseAmp'); p.backgroundNoiseAmplitude = p.backNoiseAmp; end

if isfield(e, 'sSrtTimes')
    e.stimPeriodOnOffTimes = zeros(1,length(e.sSrtTimes)+length(e.endTrialTimes)); 
    e.stimPeriodOnOffTimes(1:2:end) = e.sSrtTimes; 
    e.stimPeriodOnOffTimes(2:2:end) = e.feedbackTimes; 
    
    e.stimPeriodOnOffValues = zeros(1,length(e.sSrtTimes)+length(e.endTrialTimes)); 
    e.stimPeriodOnOffValues(1:2:end) = 1; 
end
if isfield(e, 'intOTimes')
    e.closedLoopOnOffTimes = zeros(1,length(e.intOTimes)+length(e.feedbackTimes)); 
    e.closedLoopOnOffTimes(1:2:end) = e.intOTimes; 
    e.closedLoopOnOffTimes(2:2:end) = e.feedbackTimes; 
    
    e.closedLoopOnOffValues = zeros(1,length(e.intOTimes)+length(e.feedbackTimes)); 
    e.closedLoopOnOffValues(1:2:end) = 1; 
end
if isfield(e, 'stmVValues')
    e.visStimOnOffTimes = e.stmVTimes; e.visStimOnOffValues = e.stmVValues;
    e.audStimOnOffTimes = e.stmATimes; e.audStimOnOffValues = e.stmAValues;
end

%Modify repeat on incorrect parameter field to deal with historical issues.
if ~isfield(e, 'repeatNumValues')
    e.repeatNumValues = e.endTrialTimes*0+1; 
    [v(:).maxRepeatIncorrect] = deal(0);
    p.maxRepeatIncorrect = 0;
elseif isfield(v, 'maxRetryIfIncorrect')
    tDat = {v.maxRetryIfIncorrect}'; [v.maxRepeatIncorrect] = tDat{:};
    p.maxRepeatIncorrect = p.maxRetryIfIncorrect;
elseif ~isfield(v, 'maxRepeatIncorrect')   
    [v(:).maxRepeatIncorrect] = deal(9); 
    p.maxRepeatIncorrect = 0;
end


if isfield(e, 'stimStartTimes'); warning('Stopping for debug'); keyboard; end
if ~isfield(v, 'visContrast')
    if isfield(v, 'visualContrast')     
        tDat = {v.visualContrast}'; [v.visContrast] = tDat{:};
        p.visContrast = p.visualContrast;
    elseif isfield(v, 'stimulusContrast')
        tDat = {v.stimulusContrast}'; [v.visContrast] = tDat{:};
        p.visContrast = p.stimulusContrast;
    else; [v.visContrast] = deal(1); p.visContrast = 1;
    end
end
if isfield(v, 'clickAmpDurRate')
    tDat = cellfun(@(x) x(1), {v.clickAmpDurRate}, 'uni', 0); [v.audAmplitude] = tDat{:};
    tDat = cellfun(@(x) x(2), {v.clickAmpDurRate}, 'uni', 0); [v.clickDuration] = tDat{:};
    tDat = cellfun(@(x) x(3), {v.clickAmpDurRate}, 'uni', 0); [v.clickRate] = tDat{:};
    tDat = {v.vStimAltitude}'; [v.visAltitude] = tDat{:};
    tDat = {v.vStimSigma}'; [v.visSigma] = tDat{:};

    
    p.audAmplitude = p.clickAmpDurRate(1);
    p.clickDuration = p.clickAmpDurRate(2);
    p.clickRate = p.clickAmpDurRate(3);
    p.visAltitude = p.vStimAltitude(1);
    p.visSigma = p.vStimSigma;
elseif isfield(v, 'clickDurRate')
    tDat = cellfun(@(x) x(1), {v.clickDurRate}, 'uni', 0); [v.clickDuration] = tDat{:};
    tDat = cellfun(@(x) x(2), {v.clickDurRate}, 'uni', 0); [v.clickRate] = tDat{:};
    tDat = cellfun(@(x) x(1), {v.visualAltitudeSigma}, 'uni', 0); [v.visAltitude] = tDat{:};
    tDat = cellfun(@(x) x(2), {v.visualAltitudeSigma}, 'uni', 0); [v.visSigma] = tDat{:};
    tDat = {v.audioAmplitude}'; [v.audAmplitude] = tDat{:};
    
    p.audAmplitude = p.audioAmplitude;
    p.clickDuration = p.clickDurRate(1);
    p.clickRate = p.clickDurRate(2);
    p.visAltitude = p.visualAltitudeSigma(1);
    p.visSigma = p.visualAltitudeSigma(2:3);
end


if ~isfield(p, 'audVisAzimuth') && (isfield(p, 'stimulusAzimuth') && isfield(e, 'sPosValues'))
    e.audAzimuthValues = e.sPosValues; e.audAzimuthTimes = e.sPosTimes;
    e.visAzimuthValues = e.sPosValues; e.visAzimuthTimes = e.sPosTimes;
    
    tDat = mat2cell([v.stimulusAzimuth]', ones(nTri,1));  
    [v.audInitialAzimuth] = tDat{:}; [v.visInitialAzimuth] = tDat{:};
elseif isfield(p, 'audVisAzimuth') && ~isfield(e, 'iAziValues') && ~isfield(e, 'audInitialAzimuth')
    p.audInitialAzimuth = p.audVisAzimuth(1,:); 
    p.visInitialAzimuth = p.audVisAzimuth(2,:);
    tDat = num2cell([v.audVisAzimuth]');
    [v.audInitialAzimuth] = tDat{:,1}; [v.visInitialAzimuth] = tDat{:,2};
elseif isfield(p, 'audVisAzimuth') && ~isfield(e, 'audInitialAzimuth')
    p.audInitialAzimuth = p.audVisAzimuth(1,:);
    p.visInitialAzimuth = p.audVisAzimuth(2,:);
    tDat = num2cell([e.iAziValues(1:2:end)' e.iAziValues(2:2:end)']);
    [v.audInitialAzimuth] = tDat{:,1}; [v.visInitialAzimuth] = tDat{:,2};
end

if isfield(e, 'corRValues'); tDat = num2cell(e.corRValues'); [v.correctResponse] = tDat{:};
elseif isfield(e, 'correctResponseValues'); tDat = num2cell(e.correctResponse'); [v.correctResponse] = tDat{:};
end

totalTimeOffset = b.experimentStartedTime-e.expStartTimes;
if isfield(b, 'blockTimeOffset'); totalTimeOffset = totalTimeOffset-b.blockTimeOffset(1); end

fieldList = fieldnames(e);
fieldList = fieldList(contains(fieldList, 'Times'));
for i = 1:length(fieldList); eval(['e.' fieldList{i} ' = e.' fieldList{i} ' + totalTimeOffset;']); end

fieldList = fieldnames(b.inputs);
fieldList = fieldList(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), fieldList));
for i = 1:length(fieldList); eval(['b.inputs.' fieldList{i} ' = b.inputs.' fieldList{i} ' + totalTimeOffset;']); end

if isfield(p, 'interactPunishDelays')
    p.openLoopDuration = p.interactPunishDelays(1);
    p.delayAfterIncorrect = p.interactPunishDelays(2);
else; warning('DEBUG'); keyboard;
end

if isfield(p, 'interactSigOnDurAmp')
    p.closedLoopOnsetToneAmplitude = p.interactSigOnDurAmp(3);
else; warning('DEBUG'); keyboard;
end

if isfield(p, 'rewardDurSize')
    p.delayAfterCorrect = p.rewardDurSize(1);
    p.rewardSize = p.rewardDurSize(2);
else; warning('DEBUG'); keyboard;
end

if isfield(p, 'noiseBurstAmpDur')
    p.noiseBurstAmplitude = p.noiseBurstAmpDur(1);
    p.noiseBurstDuration = p.noiseBurstAmpDur(2);
else; warning('DEBUG'); keyboard;
end

if isfield(p, 'stimulusDurRep')
    p.stimDuration = p.stimulusDurRep(1);
    p.stimContinuous = p.stimulusDurRep(2);
else; warning('DEBUG'); keyboard;
end

if isfield(p, 'preStimQuiRangeThr')
    p.preStimQuiescentRange = sort(p.preStimQuiRangeThr(1:2));
    p.preStimQuiescentThreshold = p.preStimQuiRangeThr(3);   
    tDat = num2cell(mean(p.preStimQuiescentRange)*ones(1,length(e.newTrialTimes)))';  
    [v.preStimQuiescentDuration] = tDat{:};  
elseif isfield(e, 'preStimQuiescentDurationValues')
    tDat = num2cell(e.preStimQuiescentDurationValues)'; [v.preStimQuiescentDuration] = tDat{:};
end

p = chkThenRemoveFields(p, f2Re);
v = chkThenRemoveFields(v, f2Re(~contains(f2Re, 'correctResponse')));
e = chkThenRemoveFields(e, f2Re);
b.paramsTimes = b.paramsTimes+totalTimeOffset;
b.events = e;
b.paramsValues = v;
end


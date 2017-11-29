function [standardizedBlock, standardizedParams] = standardBlkNames(block, params)
% All trials are valid trials if there are no repeats.
b = block;
e = block.events;
v = block.paramsValues;
p = params;

f2Re = {'audDevIdx';'audSampleRate';'numAudChannels';'type'; 'services'; 'defFunction'; 'experimentIdx'; 'correctResponse';...
    'servicesDescription'; 'clickAmpDurRate'; 'vStimAltitude'; 'stimulusAzimuth'; 'audVisAzimuth'; 'vStimSigma'; ...
    'interactPunishDelays'; 'stimulusDurRep'; 'noiseBurstAmpDur'; 'rewardDurSize'; 'interactSigOnDurAmp'; 'audVisThreshold'; ...
    'stimulusContrast'; 'maxRetryIfIncorrect'; 'backNoiseAmp'; 'preStimQuiRangeThr'; 'interTrialDelay'; ...
    'audioAmplitude'; 'clickDurRate'; 'visualAltitudeSigma'; 'visualContrast'; 'reflectAzimuthAndCorr'; ...
    ...
    ...
    'sPosTimes'; 'sPosValues'; 'stimStartTimes'; 'stimStartValues'; 'sSrtTimes'; 'sSrtValues'; 'stmVValues'; 'stmVTimes'; ...
    'stmAValues'; 'stmATimes'; 'intOValues'; 'intOTimes'; 'fBckValues'; 'fBckTimes'; 'rTotValues'; 'rTotTimes'; 'positiveFeedbackDuration'; ...
    'visInitialAzimuthValues'; 'visInitialAzimuthTimes'; 'audInitialAzimuthValues'; 'audInitialAzimuthTimes'; 'rewardProbabilityOnCorrect'; ...
    'preStimQuiescentDurationValues'; 'preStimQuiescentDurationTimes'; 'aPosValues'; 'aPosTimes'; 'vPosValues'; 'vPosTimes';...
    'iAziTimes'; 'iAziValues'; 'aViCValues'; 'aViCTimes'; 'visCValues'; 'visCTimes'; 'audCValues'; 'audCTimes'; 'reflectAzimuthAndCorrectResponse'; ...
    'aViMTimes'; 'aViMValues'; 'corRValues'; 'corRTimes'; 'sPreTimes'; 'sPreValues'; 'stimContinuous'; 'galvoCoordID'};

if isfield(e, 'fBckTimes'); e.feedbackTimes = e.fBckTimes; e.feedbackValues = e.fBckValues; end
if isfield(e, 'stimStartTimes'); e.sSrtTimes = e.stimStartTimes; end
if isfield(e, 'stimStartTimes'); e.sSrtTimes = e.stimStartTimes; end
if isfield(p, 'backNoiseAmp'); p.backgroundNoiseAmplitude = p.backNoiseAmp; end

if isfield(e, 'sPreValues')
    e.stimPeriodOnOffTimes = e.sPreTimes;
    e.stimPeriodOnOffValues = e.sPreValues;
elseif isfield(e, 'sSrtTimes')
    if isfield(e, 'feedbackTimes'); tDat = e.feedbackTimes; else; tDat = e.endTrialTimes; end
    e.stimPeriodOnOffTimes = zeros(1,length(e.sSrtTimes)+length(tDat));
    e.stimPeriodOnOffTimes(1:2:end) = e.sSrtTimes;
    e.stimPeriodOnOffTimes(2:2:end) = tDat;
    
    e.stimPeriodOnOffValues = zeros(1,length(e.sSrtTimes)+length(tDat));
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

    
    p.audAmplitude = p.clickAmpDurRate(1,:);
    p.clickDuration = p.clickAmpDurRate(2,:);
    p.clickRate = p.clickAmpDurRate(3,:);
    p.visAltitude = p.vStimAltitude(1,:);
    p.visSigma = p.vStimSigma;
elseif isfield(v, 'clickDurRate')
    tDat = cellfun(@(x) x(1), {v.clickDurRate}, 'uni', 0); [v.clickDuration] = tDat{:};
    tDat = cellfun(@(x) x(2), {v.clickDurRate}, 'uni', 0); [v.clickRate] = tDat{:};
    tDat = cellfun(@(x) x(1), {v.visualAltitudeSigma}, 'uni', 0); [v.visAltitude] = tDat{:};
    tDat = cellfun(@(x) x(2:3), {v.visualAltitudeSigma}, 'uni', 0); [v.visSigma] = tDat{:};
    tDat = {v.audioAmplitude}'; [v.audAmplitude] = tDat{:};
    
    p.audAmplitude = p.audioAmplitude;
    p.clickDuration = p.clickDurRate(1);
    p.clickRate = p.clickDurRate(2);
    p.visAltitude = p.visualAltitudeSigma(1);
    p.visSigma = p.visualAltitudeSigma(2:3);
end

if isfield(e, 'aPosValues')
    e.audAzimuthValues = e.aPosValues; e.audAzimuthTimes = e.aPosTimes;
    e.visAzimuthValues = e.aPosValues; e.visAzimuthTimes = e.aPosTimes;
end


if ~isfield(p, 'audVisAzimuth') && (isfield(p, 'stimulusAzimuth') && isfield(e, 'sPosValues'))
    e.audAzimuthValues = e.sPosValues; e.audAzimuthTimes = e.sPosTimes;
    e.visAzimuthValues = e.sPosValues; e.visAzimuthTimes = e.sPosTimes;
    
    tDat = mat2cell([v.stimulusAzimuth]', ones(length(e.newTrialTimes),1));  
    [v.audInitialAzimuth] = tDat{:}; [v.visInitialAzimuth] = tDat{:};
    p.audInitialAzimuth = p.stimulusAzimuth; p.visInitialAzimuth = p.stimulusAzimuth;
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
elseif isfield(e, 'preStimQuiescentDurationValues') && isfield(p, 'reflectAzimuthAndCorrectResponse')
    if  isempty(strfind(b.expDef, 'multiTemporalWorld'))
        tDat = num2cell([e.audInitialAzimuthValues' e.visInitialAzimuthValues']);
        [v.audInitialAzimuth] = tDat{:,1}; [v.visInitialAzimuth] = tDat{:,2};
    else
        tDat = num2cell(e.audInitialAzimuthValues'); [v.audInitialAzimuth] = tDat{:,1};
        tDat = num2cell(e.visContrastValues',2); [v.visContrast] = tDat{:,1};
    end
end

if isfield(e, 'corRValues'); tDat = num2cell(e.corRValues'); [v.correctResponse] = tDat{:};
elseif isfield(e, 'correctResponseValues'); tDat = num2cell(e.correctResponseValues'); [v.correctResponse] = tDat{:};
end

totalTimeOffset = b.experimentStartedTime-e.expStartTimes;
if isfield(b, 'blockTimeOffset'); totalTimeOffset = totalTimeOffset-b.blockTimeOffset(1); end

fieldList = fieldnames(e);
for i = 2:2:length(fieldList); eval(['e.' fieldList{i} ' = e.' fieldList{i} ' + totalTimeOffset;']); end

fieldList = fieldnames(b.inputs);
fieldList = fieldList(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), fieldList));
for i = 1:length(fieldList); eval(['b.inputs.' fieldList{i} ' = b.inputs.' fieldList{i} ' + totalTimeOffset;']); end

if isfield(p, 'interactPunishDelays')
    p.openLoopDuration = p.interactPunishDelays(1);
    p.delayAfterIncorrect = p.interactPunishDelays(2);
    if length(p.interactPunishDelays) > 2; p.laserDuration = p.interactPunishDelays(3); end
elseif ~isfield(p, 'laserDuration'); warning('DEBUG'); keyboard;
end

if isfield(p, 'interactSigOnDurAmp')
    p.closedLoopOnsetToneAmplitude = p.interactSigOnDurAmp(3);
elseif ~isfield(p, 'closedLoopOnsetToneAmplitude'); warning('DEBUG'); keyboard;
end

if isfield(p, 'rewardDurSize')
    p.delayAfterCorrect = p.rewardDurSize(1);
    p.rewardSize = p.rewardDurSize(2);
elseif ~isfield(p, 'rewardSize'); warning('DEBUG'); keyboard;
end

if isfield(p, 'noiseBurstAmpDur')
    p.noiseBurstAmplitude = p.noiseBurstAmpDur(1);
    p.noiseBurstDuration = p.noiseBurstAmpDur(2);
elseif ~isfield(p, 'noiseBurstDuration'); warning('DEBUG'); keyboard;
end

if isfield(p, 'stimulusDurRep')
    p.stimDuration = p.stimulusDurRep(1);
    p.stimContinuous = p.stimulusDurRep(2);
elseif ~isfield(p, 'stimContinuous'); warning('DEBUG'); keyboard;
end

if isfield(p, 'preStimQuiRangeThr')
    p.preStimQuiescentRange = sort(p.preStimQuiRangeThr(1:2));
    p.preStimQuiescentThreshold = p.preStimQuiRangeThr(3);   
    tDat = num2cell(mean(p.preStimQuiescentRange)*ones(1,length(e.newTrialTimes)))';  
    [v.preStimQuiescentDuration] = tDat{:};  
elseif isfield(e, 'preStimQuiescentDurationValues')
    tDat = num2cell(e.preStimQuiescentDurationValues)';
    if length(tDat) == length(v)-1; tDat = [tDat;0]; end
    [v.preStimQuiescentDuration] = tDat{:};
end

if ~isfield(e, 'galvoPosValues') || ~isstruct(b.galvoLog)
    b.galvoLog = struct;
    e.galvoTypeValues = 0*e.newTrialValues+1;
    e.laserTypeValues = 0*e.newTrialValues;
    e.laserPowerValues = 0*e.newTrialValues;
    e.galvoAndLaserEndTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1)+1.5;
    e.galvoCoordsValues = [0 0];
    e.galvoPosValues = 0*e.newTrialValues+1;
    p.galvoType = 1;
    p.laserPower = 0;
    p.laserTypeProportions = [1 0 0]';
    b.galvoLog.trialNum = 1:length(e.newTrialTimes);
    e.laserInitialisationTimes = b.galvoLog.trialNum*0;
    p.laserDuration = 0;
    e.galvoTTLTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1);
    [v.laserDuration] = deal(0);
end
if ~isfield(e, 'galvoTTLTimes')
    e.galvoTTLTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1);
end
if ~isfield(e, 'galvoAndLaserEndTimes')
    e.galvoAndLaserEndTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1)+e.laserDurationValues(1:sum(e.stimPeriodOnOffValues==1));
    p.laserDuration = 1.5;
end
if ~isfield(b.galvoLog, 'tictoc')
    e.laserInitialisationTimes = deal(0);
else
    if any(isnan(b.galvoLog.delay_issueLaser(b.galvoLog.laserType>0))); keyboard; end
    e.laserInitialisationTimes = b.galvoLog.tictoc;
end
if ~isfield(p, 'laserDuration'); p.laserDuration = 1.5; end
e.laserTypeValues(~ismember(1:length(e.newTrialTimes), b.galvoLog.trialNum'))=0;
p.rewardTotal = sum(p.rewardSize*e.feedbackValues>0);

paramFields = fields(p);
validConditions = p.numRepeats~=0;
for i = 1:numel(paramFields)
    if strcmp(paramFields{i}, 'type'); continue; end
    if size(p.(paramFields{i}),2) > 1
         p.(paramFields{i}) = p.(paramFields{i})(:,validConditions);
    end
end

if isfield(p, 'reflectAzimuthAndCorr'); p.reflectAzimuthAndCorrectResponse = p.reflectAzimuthAndCorr; end
if isfield(p, 'audioAmplitude'); p.audAmplitude = p.audioAmplitude; end
if isfield(p, 'reflectAzimuthAndCorrectResponse') && p.reflectAzimuthAndCorrectResponse == 1
    if isempty(strfind(b.expDef, 'multiTemporalWorld'))
        flippedIdx = max([p.visContrast;p.audAmplitude.*abs(p.audInitialAzimuth)],[],1)>0;
    else
        flippedIdx = p.coherence~=0.5;
    end
    if isfield(p, 'visInitialAzimuth') && length(p.visInitialAzimuth)==1; p.visInitialAzimuth = p.numRepeats*0+p.visInitialAzimuth; end
    if isfield(p, 'audInitialAzimuth') && length(p.audInitialAzimuth)==1; p.audInitialAzimuth = p.numRepeats*0+p.audInitialAzimuth; end
    p.numRepeats = [p.numRepeats, p.numRepeats(flippedIdx)];
end

if p.stimContinuous == 1; p.stimDuration = inf; end

for i = 1:numel(paramFields)
    if strcmp(paramFields{i}, 'type'); continue; end
        if isfield(p, 'reflectAzimuthAndCorrectResponse') && p.reflectAzimuthAndCorrectResponse == 1
            if contains(paramFields{i}, 'InitialAzimuth') || contains(paramFields{i}, 'correctResponse')
                p.(paramFields{i}) = [p.(paramFields{i}) p.(paramFields{i})(:,flippedIdx)*-1];
            elseif contains(paramFields{i}, 'coherence')
                p.(paramFields{i}) = [p.(paramFields{i}) 1-p.(paramFields{i})(:,flippedIdx)];
            elseif ~strcmp(paramFields{i}, 'numRepeats') && size(p.(paramFields{i}),2) > 1
                p.(paramFields{i}) = [p.(paramFields{i}) p.(paramFields{i})(:,flippedIdx)];
            elseif ~isempty(strfind(b.expDef, 'multiTemporalWorld')) && contains(paramFields{i}, 'visContrast')
                p.(paramFields{i}) = [p.(paramFields{i}) flip(p.(paramFields{i})(:,flippedIdx),2)];
            end
            if size(p.(paramFields{i}), 2) > 1
                p.(paramFields{i}) = p.(paramFields{i})(:,p.numRepeats>0);
                if size(unique(p.(paramFields{i})', 'rows'), 1) == 1 && ~strcmp(paramFields{i}, 'numRepeats')
                    p.(paramFields{i}) = unique(p.(paramFields{i})', 'rows')';
                end
            end
        end
end
if ~isfield(e, 'rewardAvailableValues')
    e.rewardAvailableValues = 0*e.newTrialValues+1;    
end

if sum(sum(p.laserTypeProportions(2:3,:))) == 0
    e.galvoTypeValues = 0*e.newTrialValues+1;
    e.laserTypeValues = 0*e.newTrialValues;
    e.laserPowerValues = 0*e.newTrialValues;
    e.galvoCoordsValues = [0 0];
    e.galvoPosValues = 0*e.newTrialValues+1;
    p.laserPower = 0;
    p.laserDuration = 0;
    e.laserInitialisationTimes = b.galvoLog.trialNum*0;
    e.galvoTTLTimes = e.stimPeriodOnOffTimes(e.stimPeriodOnOffValues==1);
    [v.laserDuration] = deal(0);
end

standardizedParams = prc.chkThenRemoveFields(p, f2Re);
standardizedBlock = b;
standardizedBlock.paramsValues = prc.chkThenRemoveFields(v, f2Re(~contains(f2Re, 'correctResponse')));
standardizedBlock.events = prc.chkThenRemoveFields(e, f2Re);
standardizedBlock.paramsTimes = b.paramsTimes+totalTimeOffset;
end


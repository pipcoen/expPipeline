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
    'rewardDurSize'; 'interactSigOnDurAmp'; 'audVisThreshold'};
if isfield(e, 'stimStartTimes'); e.sSrtTimes = e.stimStartTimes; end

if ~isfield(e, 'repeatNumValues')
    e.repeatNumValues = e.endTrialTimes*0+1; 
    [v(:).maxRetryIfIncorrect] = deal(0);
elseif ~isfield(v, 'maxRetryIfIncorrect')
    [v(:).maxRetryIfIncorrect] = deal(9); 
end

if isfield(e, 'stimStartTimes'); e.sSrtTimes = e.stimStartTimes; end
if ~isfield(v, 'visualContrast')
    if ~isfield(v, 'stimulusContrast')
        [v.visualContrast] = deal(1); p.visualContrast = 1;
    else; tDat = {v.stimulusContrast}'; [v.visualContrast] = tDat{:}; 
        p.visualContrast = p.stimulusContrast; p = rmfield(p, 'stimulusContrast');
    end
end
if isfield(v, 'clickAmpDurRate')
    tDat = cellfun(@(x) x(2:3), {v.clickAmpDurRate}, 'uni', 0); [v.clickDurRate] = tDat{:};
    tDat = cellfun(@(x) x(1), {v.clickAmpDurRate}, 'uni', 0); [v.audioAmplitude] = tDat{:};
    tDat = cellfun(@(x,y) [x;y], {v.vStimAltitude}, {v.vStimSigma}, 'uni', 0); [v.visualAltitudeSigma] = tDat{:};
    
    p.clickDurRate = p.clickAmpDurRate(2:3);
    p.audioAmplitude = p.clickAmpDurRate(1);
    p.visualAltitudeSigma = [p.vStimAltitude; p.vStimSigma];
end
if ~isfield(e, 'audVisAzimuth') && (isfield(p, 'stimulusAzimuth') && isfield(e, 'sPosValues'))
    e.aPosValues = e.sPosValues; e.vPosValues = e.sPosValues;
    e.aPosTimes = e.sPosTimes;  e.vPosTimes = e.sPosTimes;
    e = rmfield(e, {'sPosTimes'; 'sPosValues'});
    
    p.audioAzimuth = p.stimulusAzimuth; 
    p.visualAzimuth = p.stimulusAzimuth;
    tDat = mat2cell([v.stimulusAzimuth]', ones(nTri,1));  
    [v.audioAzimuth] = tDat{:}; [v.visualAzimuth] = tDat{:};
elseif isfield(p, 'audVisAzimuth') && ~isfield(e, 'iAziValues')
    p.audioAzimuth = p.audVisAzimuth(1,:); 
    p.visualAzimuth = p.audVisAzimuth(2,:);
    tDat = num2cell([v.audVisAzimuth]');
    [v.audioAzimuth] = tDat{:,1}; [v.visualAzimuth] = tDat{:,2};
elseif isfield(p, 'audVisAzimuth')
    p.audioAzimuth = p.audVisAzimuth(1,:);
    p.visualAzimuth = p.audVisAzimuth(2,:);
    tDat = num2cell([e.iAziValues(1:2:end)' e.iAziValues(2:2:end)']);
    [v.audioAzimuth] = tDat{:,1}; [v.visualAzimuth] = tDat{:,2};
end

if isfield(e, 'corRValues')
    tDat = num2cell(e.corRValues');
    [v.correctResponse] = tDat{:};
end

rTim = b.experimentStartedTime-e.expStartTimes;
if isfield(b, 'be2t'); rTim = rTim -b.be2t(1); end

flds = fieldnames(e);
flds = flds(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), flds));
for i = 1:length(flds); eval(['e.' flds{i} ' = e.' flds{i} ' + rTim;']); end

flds = fieldnames(b.inputs);
flds = flds(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), flds));
for i = 1:length(flds); eval(['b.inputs.' flds{i} ' = b.inputs.' flds{i} ' + rTim;']); end

if isfield(p, 'interactPunishDelays')
    p.interactDelay = p.interactPunishDelays(1);
    p.punishDelay = p.interactPunishDelays(2);
end

if isfield(p, 'interactSigOnDurAmp')
    p.interactOn = p.interactSigOnDurAmp(1);
    p.interactDuration = p.interactSigOnDurAmp(2);
    p.interactAmplitude = p.interactSigOnDurAmp(3);
end

if isfield(p, 'rewardDurSize')
    p.rewardDuration = p.rewardDurSize(1);
    p.rewardSize = p.rewardDurSize(2);
end

if isfield(p, 'noiseBurstAmpDur')
    p.noiseBurstAmplitude = p.noiseBurstAmpDur(1);
    p.noiseBurstDuration = p.noiseBurstAmpDur(2);
end

if isfield(p, 'stimulusDurRep')
    p.stimDuration = p.stimulusDurRep(1);
    p.stimContinuous = p.stimulusDurRep(2);
end

p = chkThenRemoveFields(p, f2Re);
b.paramsTimes = b.paramsTimes+rTim;
b.events = e;
b.paramsValues = v;
end


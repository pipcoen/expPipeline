function [ o ] = extractTimelineInfo( x )
%extractTimelineInfo: gets the actual differences between the time of the
%TTL pulse and the time of stimulus onset, and reports any inactivation
%trials on which the difference between the intended delay and the actual
%delay is too extreme.

t = x.timeline;
b = x.standardizedBlock;
o = struct;

name = {t.hw.inputs.name};
column = [t.hw.inputs.arrayColumn];

t.timebase = t.rawDAQTimestamps';
t.photodiode = t.rawDAQData(:,column(strcmp(name,'photoDiode')));

numTrials = length(x.validTrials);

b.trialStartTimes = b.events.trialNumTimes(1:numTrials)';
stimOnTimes = b.events.stimPeriodOnOffTimes(b.events.stimPeriodOnOffValues)';
b.stimOnTimes = stimOnTimes(1:numTrials);

t.TTL = t.rawDAQData(:,column(strcmp(name,'TTL')));
b.TTL_onsets = b.outputs.digitalTTLTimes';
b.TTL_onsets = b.TTL_onsets(1:2:end);
b.TTL_onsets = b.TTL_onsets(1:numTrials);

b.OnsetDelay = (b.TTL_onsets - b.stimOnTimes);
b.IntendedOnsetDelay = b.events.laserOnsetDelayValues(1:numTrials)';
%           figure; plot(b.IntendedOnsetDelay,b.OnsetDelay - b.IntendedOnsetDelay,'ro'); %Plot of the error at different intended delays

t.TTL = t.TTL-min(t.TTL);
t.TTL = t.TTL/max(t.TTL);
t.TTL = round(t.TTL);

t.TTL_onsets = t.timebase(logical([diff(t.TTL)==1; 0]));
t.TTL_offsets = t.timebase(logical([diff(t.TTL)==-1; 0]));
%             pulse_duration = t.TTL_offsets - t.TTL_onsets;
%             pulse_duration = round(pulse_duration,3);
%
%             t.laserOnTimes = t.TTL_onsets + 5/1000;

t.TTL_onsets = t.TTL_onsets(1:numTrials);
blockTimelineDelay = median(t.TTL_onsets - b.TTL_onsets);
o.trialStartTimes = b.trialStartTimes + blockTimelineDelay; %Get corrected trial start times
o.expectedScreenOnTime = b.stimOnTimes + blockTimelineDelay; %Get expected time of stimulus onset
o.expectedTTLTime = b.TTL_onsets + blockTimelineDelay; %Get expected time of TTL onset

%Now go through each trial, and identify the time of the laser
%and stimulus onset

photodiodeThreshold = 0.15;
ttlThreshold = 0.1;

stimOnsetPhotodiode = nan(numTrials, 1);
ttlTimeByTrial = nan(numTrials, 1);

for trial = 1:numTrials
    
    %Get photodiode and ttl trace near time when we expect the stimulus to
    %appear
    trialIndex = o.trialStartTimes(trial,1) < t.timebase & t.timebase < o.trialStartTimes(trial,1)+max(diff(o.trialStartTimes));
    
    %1) Get time around screen update
    time = t.timebase(trialIndex);
    photodiode = t.photodiode(trialIndex);
    photodiode = photodiode - photodiode(1);
    photodiode = abs(photodiode);
    photodiode = photodiode/max(photodiode);
    stimOnsetPhotodiode(trial,1) = time(find(photodiode > photodiodeThreshold,1,'first'));
    %Was originally time(ix), used to plot pd trace, may be useful
    %to know if I ever want to go back to the other code to do that
    
    %2) Get time around TTL
    ttl = t.TTL(trialIndex);
    ttlTimeByTrial(trial,1) = time(find(ttl > ttlThreshold,1,'first'));
    %Same as above with the TTL trace
end


%Peter originally had it plotting several things refer back to his version
%for those


%Calculate the average discrepency between expected stimulus time,
%and actual stimulus time. Replace every 'actual' time with
%whatever is calculated from the median. This is required because
%the photodiode measure seems unreliable sometimes
%         t.screenOnTime = t.expScreenOnTime + median( pd_stimTime - t.expScreenOnTime );
%         t.laserOnTime = t.expTTLTime + median( ttl_stimTime - t.expTTLTime ) + 5/1000;
o.screenOnTime = stimOnsetPhotodiode;
o.laserOnTime = ttlTimeByTrial + 5/1000; %Why the addition?
o.laserOnsetDelay = o.laserOnTime - o.screenOnTime;

%Remove trials where the discrepancy between real and intended onset delay times (on laser trials) is huge
try
    o.laserTrials = find(b.events.laserTypeValues>0);
    error = o.laserOnsetDelay(o.laserTrials) - b.IntendedOnsetDelay(o.laserTrials);
    o.badTrials = o.laserTrials(abs(error) > 0.1);
    o.goodTrials = setdiff(1:numTrials,o.badTrials);
%     D = structfun(@(f) f(goodTrials,:), D, 'uni', 0); %?
    fprintf('%d bad trials out of %d laser trials\n',length(o.badTrials),length(o.laserTrials));
catch me
    disp(me);
end

end


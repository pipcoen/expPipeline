function sortedByTrial = indexByTrial(trialTimes, prmTimes, prmValues, timesToSubtract, paramsForSubtraction)
%% A helper function to split variables into cells that contain all the values for that trial
%Inputs(default)
%trialTimes               (required) is an nx2 vector of the form [trialStartTimes, trialEndTimes];
%prmTimes                 (required) is an nx1 vector of times for the values that are to be split by trial.
%prmValues                (prmTimes) is an nxm vector of values that should be split.
%timesToSubtract          (0*nx1) is a nx1 of times to subtract from the paramters (e.g. if you wanted timings relative to the stimulus onset)
%paramsForSubtraction     (0*1xm) is a 1xm logical that indicates which parameter columns should have timesToSubtract subtracted

%Outputs
%sortedByTrial is a cell array of the prmValues, sorted by the trialTimes they occured between (each cell is one trial) 

%%
%Set default values
if ~exist('prmValues', 'var'); prmValues = prmTimes; end
if ~exist('timesToSubtract', 'var'); timesToSubtract =0*prmValues(1,:); end
if ~exist('paramsForSubtraction', 'var'); timesToSubtract = 0*prmTimes; end

%Use histcounts to find all the times that fally between trial start and end times--this is a computationally fast way to do this. We remove indices
%with 0 values because these are out of bounds. 
[eventCount, ~, trialIdx] = histcounts(prmTimes, sort([trialTimes(:,1);trialTimes(:,2)+realmin]));
prmValues(trialIdx==0,:) = []; trialIdx(trialIdx==0) = [];

%idxBounds finds the bounds where the trialIdx starts and ends, then removes all even indices are these lie between trials (afert end and before
%start. We also remove the eventcounts corresponding to these inter-trial spaces
idxBounds = [find(diff([-10;trialIdx])>0) find(diff([trialIdx;1e6])>0)];
idxBounds(mod(unique(trialIdx),2)==0,:) = [];
eventCount(2:2:end) = [];

%Get the unique tiral indices that aren't inter-trial spaces and pre-populate the sortedByTrial cell array according to this.
uniqueIdx = unique(trialIdx(mod(trialIdx,2)>0));
sortedByTrial = cell(length(uniqueIdx),1);
for i = 1:length(sortedByTrial)
    subtractValues = repmat(paramsForSubtraction, eventCount(i), 1);
    if any(paramsForSubtraction); subtractValues = timesToSubtract(i)*subtractValues; end

    sortedByTrial{i} = prmValues(idxBounds(i,1):idxBounds(i,2),:);
    if isempty(sortedByTrial{i}); continue; end
    sortedByTrial{i} = sortedByTrial{i} - subtractValues;
end
end
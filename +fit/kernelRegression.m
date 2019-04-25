function [fitKernels, chosenClusters, predictedSignals, cvErr] = kernelRegression(blk, events, kernalOpt)
%%
% -- spikeTrains is nS by nTimePoints, any number of signals to be fit
% -- t is 1 by nTimePoints
% -- eventTimes is a cell array of times of each event
% -- eventValues is a cell array of "values" for each event, like if you want
% different instances of the event to be fit with a scaled version of the
% kernel. E.g. contrast of stimulus or velocity of wheel movement.
% Note that this is different from a second stage estimation of cofficents
% that one might want to fit on the data, by keeping the kernels fixed.
% -- kernalOpt.windows is a cell array of 2 by 1 kernalOpt.windows, [startOffset endOfsampleRateet]
% -- kernalOpt.lambda is a scalar, the regularization amount. set to 0 to do no
% regularization (if non-zero, does ridge regression)
% -- kernalOpt.cvFold is 2 by 1, [foldSize, nToCalculate]. So [5 5] does 5-fold CV
% and calculates the error on all five of the test sets. [5 1] still holds
% out 20% but only does this once. 
% This is mainly inspired by the paper from Park and Pillow, NatureNeuro, 2013
% Armin Lak, 2017, based on initial codes from Kenneth and Nick.

if ~exist('kernalOpt', 'var'); kernalOpt = struct; end
if ~isfield(kernalOpt, 'lambda'); kernalOpt.lambda = 0; end
if ~isfield(kernalOpt, 'cvFold'); kernalOpt.cvFold = [0 0]; end
if ~isfield(kernalOpt, 'validTrialsOnly'); kernalOpt.validTrialsOnly = 1; end
if ~isfield(kernalOpt, 'eventValues'); kernalOpt.eventValues = repmat({[]}, 1, length(events)); end
if ~isfield(kernalOpt, 'windows'); kernalOpt.windows = repmat({[-0.1;0.5]}, 1, length(events)); end
if ~isfield(kernalOpt, 'zscore'); kernalOpt.zscore = 1; end
if ~isfield(kernalOpt, 'minSpikes'); kernalOpt.minSpikes = 200; end

%%
eventTimes = cell(length(events),1);
for i = 1:length(events)
    if strcmp(events{i}(end-1:end), 'On') && isfield(blk, [events{i} 'Off'])
        eventTimes{i} = blk.([events{i} 'Off']);
        if iscell(eventTimes{i}); eventTimes{i} = cell2mat(eventTimes{i}); end
        eventTimes{i} = eventTimes{i}(:,1);
    elseif strcmp(events{i}(end-2:end), 'Off') && isfield(blk, [events{i}(1:end-2) 'OnOff'])
        eventTimes{i} = blk.([events{i}(1:end-2) 'OnOff']);
        if iscell(eventTimes{i}); eventTimes{i} = cell2mat(eventTimes{i}); end
        eventTimes{i} = eventTimes{i}(:,2);
    else, eventTimes{i} = blk.(events{i});
        if iscell(eventTimes{i}); eventTimes{i} = cell2mat(eventTimes{i}); end
    end
end

%%
[spikeTrains, t, chosenClusters] = kil.mkSpikeTrains(double(blk.ephSpikeTimes), double(blk.ephSpikeClusters), 10/1000);
if kernalOpt.validTrialsOnly
    windows2Keep = [blk.trialStartEnd(:,1) blk.trialStartEnd(:,2)+0.5];
    validTimes = cell2mat(arrayfun(@(x,y) t(t>x & t<y)', windows2Keep(:,1), windows2Keep(:,2), 'uni', 0));
    validIdx = ismember(t, validTimes);
    spikeTrains = spikeTrains(:, validIdx);
    t = t(validIdx);
end

numOfSpikes = sum(spikeTrains,2);
chosenClusters = chosenClusters(numOfSpikes>kernalOpt.minSpikes);
spikeTrains = spikeTrains(numOfSpikes>kernalOpt.minSpikes,:);
if kernalOpt.zscore; spikeTrains = zscore(spikeTrains, [], 2); end
%%
sampleRate = 1/mode(diff(t));
numOfTimepoints = length(t);
numOfSignals = size(spikeTrains,1);

% this is the function used to evaluate the cross validated error. Should
% return numOfSignals x 1, the performance on each signal to be predicted.
% 2 possible functions (often give comparable results):
% cvEvalFunc = @(pred, actual)1-mean(mean((pred-actual).^2))/mean(mean(actual.^2));
cvEvalFunc = @(pred, actual)1- var(pred-actual); % here assume variance of actual is 1 - it is (or is close) if data were zscored. Otherwise you'd want to divide by it

startOffset = cellfun(@(x) x(1)*sampleRate,kernalOpt.windows);
numWindowSamps = cellfun(@(x) round(diff(x)*sampleRate),kernalOpt.windows);
numWinSampsTotal = sum(numWindowSamps);
cumsumWindows = cumsum([0 numWindowSamps]);

%%
A = zeros(numOfTimepoints,numWinSampsTotal);

for i = 1:length(eventTimes)
    [theseET, sortI] = sort(eventTimes{i}(:));
    timeIdx = prc.nearestPoint(theseET, t);
    eventFrames = round(timeIdx+startOffset(i)); % the "frames", i.e. indices of t, at which the start of each event's window occurs
    
     if isempty(kernalOpt.eventValues{i}); theseEventValues = ones(size(eventFrames));
     else, theseEventValues = kernalOpt.eventValues{i}(sortI);
     end
     
    for j = 1:numWindowSamps(i)
        theseSamps = eventFrames+j;
        inRange = theseSamps>0&theseSamps<=size(A,1);
        A(theseSamps(inRange),cumsumWindows(i)+j) = theseEventValues(inRange);
    end
    
    if kernalOpt.lambda>0
        spikeTrains(:,end+1:end+numWindowSamps(i)) = 0;
        A(end+1:end+numWindowSamps(i),cumsumWindows(i)+1:cumsumWindows(i)+numWindowSamps(i)) = diag(kernalOpt.lambda*ones(1,numWindowSamps(i)));
    end
end

%%
if kernalOpt.cvFold(1)>0
    
    cvp = cvpartition(numOfTimepoints,'KFold', kernalOpt.cvFold(1));
    cvErr = zeros(numOfSignals,kernalOpt.cvFold(2));
    for k = 1:kernalOpt.cvFold(2)
        fprintf(1, 'kernalOpt.cvFold %d/%d\n', k, kernalOpt.cvFold(2))
        if kernalOpt.lambda>0
            % if using regularization, you want the regularization rows to
            % always be part of training
            trainInds = vertcat(cvp.training(k), true(size(spikeTrains,2)-numOfTimepoints,1));
        else
            trainInds = cvp.training(k);
        end
        
        testInds = cvp.test(k);
        
        trainSetObservations = spikeTrains(:,trainInds);
        trainSetPredictors = A(trainInds,:);
        X = solveLinEq(trainSetPredictors,trainSetObservations'); % X becomes numWinSampsTotal by nS 
    
        predictedSignals = (A(testInds,:)*X);
        testSetObservations = spikeTrains(:,testInds)';
        cvErr(:,k) = cvEvalFunc(predictedSignals, testSetObservations);
    end
    
else
    
    X = solveLinEq(A,spikeTrains'); % X becomes numWinSampsTotal by nS
    cvErr = [];
    
end

% set outcome of empty events to Nan
for i = 1:length(eventTimes)
    if ~isempty(eventTimes{i})
    fitKernels{i} = X(cumsumWindows(i)+1:cumsumWindows(i+1),:);
    else
         fitKernels{i}=nan;  
    end
end

predictedSignals = [];
if nargout>1
    % return the predicted signal, given the kernels
    predictedSignals = (A(1:numOfTimepoints,:)*X)';
end


function X = solveLinEq(A, B)
% This is just mldivide, but it turns out to be faster, empirically, to
% make the variables gpuArrays and use pinv instead. 

%gA = gpuArray(single(A));
%gB = gpuArray(single(B));
%X = gather(pinv(gA)*gB);

% if you don't have gpu, use this one:
X = A\B;
function resp = driftDiffusion(fitOpt)
if~exist('fitOpt', 'var'); fitOpt = []; end
if~isfield(fitOpt, 'numberOfTrials'); fitOpt.numberOfTrials = 1000; end
if~isfield(fitOpt, 'upperThresh'); fitOpt.upperThresh = 1; end
if~isfield(fitOpt, 'lowerThresh'); fitOpt.lowerThresh = -1; end
if~isfield(fitOpt, 'startPnt'); fitOpt.startPnt = 0; end
if~isfield(fitOpt, 'driftRate'); fitOpt.driftRate = 0.3; end
if~isfield(fitOpt, 'stepSize'); fitOpt.stepSize = 0.01; end
if~isfield(fitOpt, 'timeLimit'); fitOpt.timeLimit = 10; end
% addParamValue(p, 'a',0.16 ); %upper threshold
% addParamValue(p, 'b', []); %lower threshold. If []. b=0
% addParamValue(p, 'z', []); %starting point, if [] then z=a/2

% addParamValue(p, 'v', 0.3); %drift rate within trial
% addParamValue(p, 'Ter', .26);  %Non decision time
% addParamValue(p, 'st', 0); %variability in the non decision time
% addParamValue(p, 'eta', 0.063);  %variability in drift rate across trial
% addParamValue(p, 'sz', 0);  %variability in starting point
% addParamValue(p, 'c', 0.1); %std within trial, put [] is you want it to calculate for you
% addParamValue(p, 'tau', 0.0001); %step

% mu=a/v; sig=(a.^2./(c^2));
resp.choice = nan(fitOpt.numberOfTrials,1);
resp.walk = cell(fitOpt.numberOfTrials,1);

for i=1:fitOpt.numberOfTrials
    randValues = normrnd(fitOpt.driftRate,0);
    timeseriesSUM=cumsum([fitOpt.startPnt; normrnd(randValues*fitOpt.stepSize,sqrt(fitOpt.stepSize),fitOpt.timeLimit/fitOpt.stepSize,1)]);
    decisionPnt = find(timeseriesSUM>fitOpt.upperThresh | timeseriesSUM<fitOpt.lowerThresh, 1);
    
    if isempty(decisionPnt)
        resp.walk{i} = timeseriesSUM;
        continue
    end
    resp.choice(i) = sign(timeseriesSUM(decisionPnt));
    resp.walk{i} = timeseriesSUM(1:decisionPnt);
end
end
%}
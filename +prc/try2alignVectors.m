function [t1Corrected, t2Corrected] = try2alignVectors(t1, t2, diffThresh)

%%
t1Diff = diff(t1);
t2Diff = diff(t2);
if ~exist('diffThresh', 'var')
    [~, dist] = knnsearch(t2Diff,t1Diff);
    diffThresh = max([0.25 20*mad(dist)]);
end

buff = 5;
offsetTest = cell2mat(arrayfun(@(x) [sum(abs(t1Diff((buff+1:buff*2)+x)-t2Diff(buff+1:buff*2))) x],-buff:buff, 'uni', 0)');
initialOffset = offsetTest(offsetTest(:,1)==min(offsetTest(:,1)),2);
if initialOffset < 0; t2(1:abs(initialOffset)) = [];
elseif initialOffset > 0; t1(1:initialOffset) = [];
end

minL = min([length(t1) length(t2)]);
compareVect = [t1(1:minL)-(t1(1)) t2(1:minL)-t2(1)];
loopNumber = 0;
%%
while find(abs(diff(diff(compareVect,[],2)))>diffThresh,1)
    errPoint = find(abs(diff(diff(compareVect,[],2)))>diffThresh,1);
    t2(errPoint+1) = [];
    t2(errPoint-1:errPoint+1) = [];
    t1(errPoint-1:errPoint+1) = [];
    
    minL = min([length(t1) length(t2)]);
    compareVect = [t1(1:minL)-(t1(1)) t2(1:minL)-t2(1)];
    loopNumber = loopNumber+1;
    if loopNumber > 50; error('Extreme alignment error'); end
end
t1Corrected = t1(1:minL);
t2Corrected = t2(1:minL);
end
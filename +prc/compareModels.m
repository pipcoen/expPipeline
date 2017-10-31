function modelErrors = compareModels
%% Function to combine and filter blk files.
s = spatialAnalysis({'PC011'; 'PC012';'PC013';'PC015'}, {{'rng', '2017-06-14','2017-08-16'}; 'last40'; 'last40'; 'last40'});
%%
for i = 1:4
    blk = s.blocks{i};
    numRightTurns = prc.makeGrid(blk, blk.response==2, @sum, 1);
    numTrials = prc.makeGrid(blk, blk.response, @length, 1);
    numLeftTurns = numTrials-numRightTurns;
    
    responseData = numRightTurns./numTrials;
    trialType = prc.makeGrid(blk, blk.trialType, @mean, 1);
    audIdx = trialType==1 | trialType==0;
    visIdx = trialType==2| trialType==0;
    
    numeratorTurns = numRightTurns(audIdx)*numRightTurns(visIdx)';
    denomenatorTurns = numLeftTurns(audIdx)*numLeftTurns(visIdx)';
    mulModel = numeratorTurns./(numeratorTurns+denomenatorTurns);
    
    modelErrors{1}(i,1) = sum(abs(responseData(2,6:8)-responseData(1,6:8))) + sum(abs(responseData(3,2:4)-responseData(2,2:4)));
    modelErrors{2}(i,1) = sum(abs(responseData(2,5)-responseData(1,6:8))) + sum(abs(responseData(3,2:4)-responseData(2,5)));
    modelErrors{3}(i,1) = sum(abs(mulModel(1,5:7)-responseData(1,6:8))) + sum(abs(mulModel(3,1:3)-responseData(3,2:4)));
end

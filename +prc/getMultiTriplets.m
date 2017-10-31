function triplets = getMultiTriplets(blk, coherence)
%% Function to combine and filter block files.
if ~exist('coherence', 'var'); coherence = 1; end
if coherence; tTyp = 3; else, tTyp = 4; end
trialTypes = prc.makeGrid(blk, blk.trialType, @mean, 3);
multiTrials = trialTypes == tTyp;
percentCorrect = prc.makeGrid(blk, blk.feedback>0, @mean, 2);
if ~coherence 
    followAuditory = prc.makeGrid(blk, (((blk.response==2)*2)-1)==sign(blk.audDiff), @mean, 2);
    percentCorrect(trialTypes == tTyp) = followAuditory(trialTypes == tTyp);
end
triplets = struct;
idx = 0;
for j = find(multiTrials)'
    vIdx = blk.grids.visValues==abs(blk.grids.visValues(j)) & trialTypes==2;
    aIdx = blk.grids.audValues==abs(blk.grids.audValues(j)) & trialTypes==1;
    if ~any(vIdx(:)) || ~any(aIdx(:)); continue; end
    triplets.yData{idx+1,1} = percentCorrect{aIdx};
    triplets.yData{idx+3,1} = percentCorrect{vIdx}*coherence + (1-(percentCorrect{vIdx}))*~coherence;
    triplets.yData{idx+2,1} = percentCorrect{j};
    triplets.xLabels(idx+1:idx+3,:) = {['Aud \newline' num2str(blk.grids.audValues(j))]; 'Mul'; ['Vis \newline' num2str((blk.grids.visValues(j)))]};
    triplets.xPosition(idx+1:idx+3,1) = floor((idx+1)/3)+(idx+1:idx+3);
    [~, pairedTestSig(1,1)] = ttest(triplets.yData{idx+1}, triplets.yData{idx+2});
    [~, pairedTestSig(2,1)] = ttest(triplets.yData{idx+2}, triplets.yData{idx+3});
    triplets.significance(floor((idx+1)/3)*2+(1:2),1) = pairedTestSig;
    triplets.testedPairs(floor((idx+1)/3)*2+(1:2),:) = {triplets.xPosition([idx+1 idx+2]); triplets.xPosition([idx+2 idx+3])};
    idx = idx+3;
end

function triplets = getCoherentConflictPairs(blk)
%% Function to combine and filter block files.
trialTypes = prc.makeGrid(blk, blk.trialType, @mean);
coherentTrials = trialTypes == 3;
reactionTimes = prc.makeGrid(blk, blk.timeToWheelMove, @mean);
numTrials = prc.makeGrid(blk, blk.responseMade==2, @length, 0);
triplets = struct;
idx = 0;
multiTrials(:,floor(size(multiTrials,2)/2):floor(size(multiTrials,2)/2)+2) = 0;
multiTrials(:,[1 end]) = 0;
multiTrials(1,:) = 0;
for j = find(multiTrials)'
    if coherence; vIdx = blk.grids.visValues==abs(blk.grids.visValues(j)) & trialTypes==2;
    else, vIdx = blk.grids.visValues==(-1*(abs(blk.grids.visValues(j)))) & trialTypes==2;
    end
    aIdx = blk.grids.audValues==abs(blk.grids.audValues(j)) & trialTypes==1;
    if ~any(vIdx(:)) || ~any(aIdx(:)); continue; end
    validDays = numTrials{aIdx}>3 & numTrials{vIdx}>3 & numTrials{j}>3;
    triplets.yData{idx+1,1} = fractionRight{aIdx}(validDays);
    triplets.yData{idx+3,1} = fractionRight{vIdx}(validDays);
    triplets.yData{idx+2,1} = fractionRight{j}(validDays);
    triplets.xTickLabels(idx+1:idx+3,:) = {['Aud \newline' num2str(blk.grids.audValues(j))]; 'Mul'; ['Vis \newline' num2str((blk.grids.visValues(j)))]};
    triplets.xTickLocations(idx+1:idx+3,1) = floor((idx+1)/3)+(idx+1:idx+3);
    triplets.pairs2test(floor((idx+1)/3)*2+(1:2),:) = {[idx+1 idx+2]; [idx+2 idx+3]};
    idx = idx+3;
end
triplets.linkedGroups = triplets.pairs2test;
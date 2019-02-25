function switchedBlock = switchBlockToContraIpsi(block)
rightIdx = block.visDiff>0 | (block.visDiff==0 & block.audDiff>0);
block.responseMade(rightIdx) = (block.responseMade(rightIdx)*-1+3).*block.responseMade(rightIdx)>0;
block.correctResponse(rightIdx) = (block.correctResponse(rightIdx)*-1+3).*block.correctResponse(rightIdx)>0;
block.audInitialAzimuth(rightIdx) = block.audInitialAzimuth(rightIdx)*-1;
block.audDiff(rightIdx) = block.audDiff(rightIdx)*-1;
block.audValues = unique(block.audDiff);

block.visInitialAzimuth(rightIdx) = block.visInitialAzimuth(rightIdx)*-1;
block.visInitialAzimuth(isinf(block.visInitialAzimuth)) = inf;
block.visDiff(rightIdx) = block.visDiff(rightIdx)*-1;
block.visValues = unique(block.visDiff);
block.galvoPosition(rightIdx,1) = -1*block.galvoPosition(rightIdx,1);


allConditions = [block.audAmplitude block.visContrast block.audInitialAzimuth block.visInitialAzimuth];
uniqueConditions = unique(allConditions, 'rows');
%Create a set of unique conditions, where each row is a condition in the order: [zero conditions; right conditions; left conditions].
leftInitialConditions = uniqueConditions(uniqueConditions(:,end-1)< 0 | ((isinf(uniqueConditions(:,end-1)) | ~(uniqueConditions(:,end-1))) & uniqueConditions(:,end)<0),:);
zeroConditions = uniqueConditions(all([any(uniqueConditions(:,[1,3])==0,2) any(uniqueConditions(:,[2,4])==0,2)],2),:);
rightInitialConditions = [leftInitialConditions(:,1:2) leftInitialConditions(:,end-1:end)*-1];
rightInitialConditions(isinf(rightInitialConditions)) = inf;
rightInitialConditions = [rightInitialConditions; uniqueConditions(~ismember(uniqueConditions, [zeroConditions; leftInitialConditions; rightInitialConditions], 'rows'),:)];
uniqueConditions = [zeroConditions; rightInitialConditions; leftInitialConditions];
[~, rightConditionsIdx] = ismember(allConditions, rightInitialConditions, 'rows');
[~, leftConditionsIdx] = ismember(allConditions, leftInitialConditions, 'rows');
if any(all([rightConditionsIdx~=0, leftConditionsIdx~=0],2)); error('Detect same condition as being Left and Right'); end
conditionLabel = rightConditionsIdx + -1*leftConditionsIdx;
uniqueConditionLabels = [0*find(zeroConditions,1); (1:size(rightInitialConditions,1))'; -1*(1:size(leftInitialConditions,1))'];
[~, conditionRowIdx] = ismember(conditionLabel, uniqueConditionLabels);
uniqueDiff = [uniqueConditions(:,3) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];

block.uniqueConditions = uniqueConditions;
block.uniqueDiff = uniqueDiff;
block.uniqueConditionLabels = uniqueConditionLabels;
block.conditionLabelRow = [conditionLabel conditionRowIdx]; 

switchedBlock = block;
end
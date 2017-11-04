function results = outerProduct(block)
numRightTurns = prc.makeGrid(block, block.responseMade==2, @sum, 1);
numTrials = prc.makeGrid(block, block.responseMade, @length, 1);
numLeftTurns = numTrials-numRightTurns;

responseData = numRightTurns./numTrials;    
results.svd.visValues = block.visValues(~any(isnan(responseData)))*100;
results.mul.visValues = block.visValues*100;

results.svd.model{1,1} = responseData(:,~any(isnan(responseData),1));
[U,S,V] = svd(results.svd.model{1,1}, 'econ');
results.svd.model{2,1} = U(:,1)*V(:,1)'*S(1,1);
results.svd.model{3,1} = results.svd.model{1,1} - results.svd.model{2,1};

audIdx = block.grids.visValues==0;
visIdx = block.grids.audValues==0;

numeratorTurns = numRightTurns(audIdx)*numRightTurns(visIdx)';
denomenatorTurns = numLeftTurns(audIdx)*numLeftTurns(visIdx)';

results.mul.model{1,1} = responseData;
results.mul.model{2,1} = numeratorTurns./(numeratorTurns+denomenatorTurns);
results.mul.model{2,1}(audIdx) = nan; results.mul.model{2,1}(visIdx) = nan;
results.mul.model{3,1} = results.mul.model{1,1} - results.mul.model{2,1};
results.mul.model{3,1}(isnan(responseData)) = nan;
end

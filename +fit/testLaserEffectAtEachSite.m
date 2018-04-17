function [deltaRightGrid, pValGrid, gridXY, nTrials] = testLaserEffectAtEachSite(block, testType, siteIdx)
switch testType
    case 'VisUni(L-R)'
        leftBlock = prc.combineBlocks(block, block.trialType==2 & block.correctResponse==1);
        rightBlock = prc.combineBlocks(block, block.trialType==2 & block.correctResponse==2);
    case 'AudUni(L-R)'
        leftBlock = prc.combineBlocks(block, block.trialType==1 & block.correctResponse==1);
        rightBlock = prc.combineBlocks(block, block.trialType==1 & block.correctResponse==2);
    case 'CohUni(VL-VR)'
        leftBlock = prc.combineBlocks(block, block.trialType==3 & block.correctResponse==1);
        rightBlock = prc.combineBlocks(block, block.trialType==3 & block.correctResponse==2);
    case 'ConUni(AL-AR)'
        leftBlock = prc.combineBlocks(block, block.trialType==4 & block.audInitialAzimuth<0);
        rightBlock = prc.combineBlocks(block, block.trialType==4 & block.audInitialAzimuth>0);
end
nShuffles = 500;
inactiveGrid = cell(nShuffles+1, 2);
for i = 1:nShuffles+1
    for j = 1:2
    if j==1; tempBlock = leftBlock; else; tempBlock = rightBlock; end
    if i > 1
        tempBlock.laserType = tempBlock.laserType(randperm(length(tempBlock.laserType))); 
        tempBlock.galvoPosition = tempBlock.galvoPosition(randperm(length(tempBlock.laserType)),:);
    end
    normBlock = prc.combineBlocks(tempBlock, tempBlock.laserType==0);
    laserBlock = prc.combineBlocks(tempBlock, tempBlock.laserType==1);
    [inactiveGrid{i,j}, gridXY] = prc.makeGrid(laserBlock, laserBlock.responseMade==2, @mean, 'galvouni');
    inactiveGrid{i,j} = inactiveGrid{i,j} - mean(normBlock.responseMade==2);
    end
end 
nTrials = prc.makeGrid(laserBlock, laserBlock.responseMade==2, @length, 'galvouni');
leftResults = reshape(cell2mat(inactiveGrid(:,1)'), [size(gridXY{1}), nShuffles+1]);
rightResults = reshape(cell2mat(inactiveGrid(:,2)'), [size(gridXY{1}), nShuffles+1]);
diffGrid = leftResults + rightResults;
%%
pValGrid = nan*diffGrid(:,:,1);
for i = 1:size(diffGrid,1)
    for j = 1:size(diffGrid,2)
        if isnan(diffGrid(i,j)); continue; end
        pValGrid(i,j) = (find(sort(abs(diffGrid(i,j,:)), 'descend')==abs(diffGrid(i,j,1)))-1)./nShuffles;
    end
end
deltaRightGrid = diffGrid(:,:,1);

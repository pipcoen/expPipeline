function separatedBlocks = separateBlockIntoInactivationSites(block)
siteList = unique([abs(block.galvoPosition(:,1)), block.galvoPosition(:,2)], 'rows');
separatedBlocks.blocks = cell(size(siteList, 1), 3);
laserType = block.laserType;
galvoPosition = block.galvoPosition;
for i = 1:size(siteList, 1)
    separatedBlocks.blocks{i,1} = combineAndFilterBlocks(block, ismember(galvoPosition, [siteList(i,1)*-1, siteList(i,2)], 'rows') & laserType==1);
    separatedBlocks.blocks{i,2} = combineAndFilterBlocks(block, ismember(galvoPosition, [siteList(i,1), siteList(i,2)], 'rows') & laserType==2);
    separatedBlocks.blocks{i,3} = combineAndFilterBlocks(block, ismember(galvoPosition, [siteList(i,1), siteList(i,2)], 'rows') & laserType==1);
end
separatedBlocks.siteList = siteList;
end
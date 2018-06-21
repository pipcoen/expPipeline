function b = getDefinedSubset(b, subsetTag)
minVis = 0.04;
switch lower(subsetTag)
    case 'vl'; b = prc.combineBlocks(b, b.trialType==2 & b.correctResponse==1 & b.visContrast>minVis);
    case 'vr'; b = prc.combineBlocks(b, b.trialType==2 & b.correctResponse==2 & b.visContrast>minVis);
    case 'al'; b = prc.combineBlocks(b, b.trialType==1 & b.audInitialAzimuth<0);
    case 'ar'; b = prc.combineBlocks(b, b.trialType==1 & b.audInitialAzimuth>0);
    case 'cohl'; b = prc.combineBlocks(b, b.trialType==3 & b.correctResponse==1 & b.visContrast>minVis);
    case 'cohr'; b = prc.combineBlocks(b, b.trialType==3 & b.correctResponse==2 & b.visContrast>minVis);
    case 'conl'; b = prc.combineBlocks(b, b.trialType==4 & b.visInitialAzimuth>0 & b.visContrast>minVis);
    case 'conr'; b = prc.combineBlocks(b, b.trialType==4 & b.visInitialAzimuth<0 & b.visContrast>minVis);
end
end
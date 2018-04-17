function b = getDefinedSubest(b, subsetTag)
switch lower(subsetTag)
    case 'vl'; b = prc.combineBlocks(b, b.trialType==2 & b.correctResponse==1 & b.visContrast > 0.08);
    case 'vr'; b = prc.combineBlocks(b, b.trialType==2 & b.correctResponse==2 & b.visContrast > 0.08);
    case 'al'; b = prc.combineBlocks(b, b.trialType==1 & b.correctResponse==1);
    case 'ar'; b = prc.combineBlocks(b, b.trialType==1 & b.correctResponse==2);
    case 'cohl'; b = prc.combineBlocks(b, b.trialType==3 & b.correctResponse==1 & b.visContrast > 0.08);
    case 'cohr'; b = prc.combineBlocks(b, b.trialType==3 & b.correctResponse==2 & b.visContrast > 0.08);
    case 'conl'; b = prc.combineBlocks(b, b.trialType==4 & b.visInitialAzimuth>0 & b.visContrast > 0.08);
    case 'conr'; b = prc.combineBlocks(b, b.trialType==4 & b.visInitialAzimuth<0 & b.visContrast > 0.08);
end

end
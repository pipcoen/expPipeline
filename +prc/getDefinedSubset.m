function b = getDefinedSubset(b, subsetTag)
minVis = 0.01;
switch lower(subsetTag)
    case 'vl'; b = prc.filtStruct(b, b.trialType==2 & b.correctResponse==1 & b.visContrast>minVis);
    case 'vr'; b = prc.filtStruct(b, b.trialType==2 & b.correctResponse==2 & b.visContrast>minVis);
    case 'al'; b = prc.filtStruct(b, b.trialType==1 & b.audInitialAzimuth<0);
    case 'ar'; b = prc.filtStruct(b, b.trialType==1 & b.audInitialAzimuth>0);
    case 'cohl'; b = prc.filtStruct(b, b.trialType==3 & b.correctResponse==1 & b.visContrast>minVis);
    case 'cohr'; b = prc.filtStruct(b, b.trialType==3 & b.correctResponse==2 & b.visContrast>minVis);
    case 'conl'; b = prc.filtStruct(b, b.trialType==4 & b.visInitialAzimuth<0 & b.visContrast>minVis);
    case 'conr'; b = prc.filtStruct(b, b.trialType==4 & b.visInitialAzimuth>0 & b.visContrast>minVis);
end
end
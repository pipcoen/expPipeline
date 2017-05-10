function blkTConv = alignBlockTimes(b, t)
diodeCha = t.hw.inputs(strcmp({t.hw.inputs.name},'photoDiode')).arrayColumn;
diodeCha = t.rawDAQData(:,diodeCha);
tLineTim = t.rawDAQTimestamps';
[~, lowHiVal] = kmeans(diodeCha, 2);
flipPnts = abs(diff(pCha > mean(lowHiVal))) > 0; % look for the flips
flipPnts = mean([pTim([tLineTim ; flipPnts]) tLineTim([flipPnts ; false])], 2);
blkFlips = b.stimWindowUpdateTimes;

if length(blkFlips)-length(flipPnts) > 10
    warning('Missed some flips on Photodiode. Will try to compensate');
    if mean(flipPnts(1:10) - blkFlips(1:10)) < 0
        error('Could not compensate');
    end
    [covarVal, covarLag] = xcov(diff(flipPnts(15:50)), diff(blkFlips(15:50)), 50);
    optimLag = covarLag(covarVal==max(covarVal));
    if optimLag < 0;  blkFlips(1:optimLag*-1)=[];
    elseif optimLag>0; flipPnts(1:optimLag) = [];
    end
    
    initDiff = mean(flipPnts(25:100) - blkFlips(25:100));
    blkFlips = blkFlips + initDiff;
    blkFlips = blkFlips(knnsearch(blkFlips, flipPnts));
    blkFlips = blkFlips - initDiff;
end

[covarVal, covarLag] = xcov(diff(flipPnts), diff(blkFlips), 50);
optimLag = covarLag(covarVal==max(covarVal));
if optimLag < 0;  blkFlips(1:optimLag*-1)=[];
elseif optimLag>0; flipPnts(1:optimLag) = [];
end

minFlips = min([length(blkFlips), length(flipPnts)]);
[blkTConv, stat] = robustfit(blkFlips(2:minFlips-2), flipPnts(2:minFlips-2));
if stat.ols_s > 0.02
    warning(['Significant error in fit: ' num2str(stat.ols_s)]);
end
end
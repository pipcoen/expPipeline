function oPut = splitByTrials(stEn, pTim, pVal, sSrt)
if nargin < 4; sSrt = 0*stEn(:,1); end
tIdx = cell2mat(fun.map(@(r) r+0*pTim(pTim>stEn(r,1) & ...
    pTim<=stEn(r,2))', 1:size(stEn,1)))';
oPut = cell(max(tIdx),1);
for i = 1:size(stEn,1)
    oPut{i,1} = pVal(tIdx==i,:) - sSrt(i);
end


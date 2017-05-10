function oPut = indexByTrial(blk, pTim, pVal, sSub)
if nargin < 4; sSub = 0*pVal(1,:); end
tIdx = cell2mat(fun.map(@(r) r+0*pTim(pTim>=blk.stEn(r,1) & ...
    pTim<=blk.stEn(r,2))', 1:size(blk.stEn,1)))';
oPut = cell(max(tIdx),1);
for i = 1:length(oPut)
    tSub = (sSub==1)*blk.sSrt(i);
    tVal = pVal(tIdx==i,:);
    if isempty(tVal); continue; end
    tSub(sSub==2) = tVal(1,sSub==2);
    oPut{i} = single(bsxfun(@minus, tVal, tSub));
end
end
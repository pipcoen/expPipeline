function spkR = normStimSpks(stEn, spkT, preT)
if ~exist('preT', 'var'); preT = 0.5; end

if ~iscell(spkT); spkT = {spkT}; end
triN = cell2mat(fun.map(@(x) countByTrials(x, stEn), spkT));
preN = cell2mat(fun.map(@(x) countByTrials(x, [stEn(:,1)-preT stEn(:,1)]), spkT));
spkR = triN - preN;
end
function spkR = countByTrials(spkT, stEn)
spkR = arrayfun(@(x) sum(spkT>stEn(x,1) & spkT<stEn(x,2))./...
    diff(stEn(x,:)), 1:size(stEn,1));
end
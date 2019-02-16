function [spikeTrains, binTimes, templates] = mkSpikeTrains(spikeTimes, spikeTemplates, windowSize)
%%
templates = unique(spikeTemplates);
binBoundaries = (min(spikeTimes):windowSize:ceil(max(spikeTimes)))';
spikeTrains = arrayfun(@(x) histc(spikeTimes(spikeTemplates == templates(x)), binBoundaries), 1:length(templates), 'uni', 0);
spikeTrains = cell2mat(spikeTrains)';
binTimes = (binBoundaries(1:end-1)+windowSize/2)';
end









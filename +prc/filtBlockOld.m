function filtered = filtBlock(unfiltered, criterion, filterTag)
%% Function to fulter structure.
if ~exist('criterion', 'var'); filtered = []; return; end
if ~exist('filterTag', 'var')
    numDat = unfiltered.number;
    numFields = fields(numDat);
    numbers = cellfun(@(x) numDat.(x), numFields(~contains(numFields, 'Per')));
    if length(numbers) ~= length(unique(numbers)); error('Please specify tag. Type of filter unclear'); end
    filterTag = numFields(numbers==length(criterion));
end
filtered = unfiltered;
criterion = criterion>0;
blockFields = fields(unfiltered);

if strcmpi(filterTag, 'experiments')
    expStructures = blockFields(contains(blockFields, {'subject', 'expDate', 'expNum', 'rigName', 'expType', 'expDef',...
        'performanceAVM', 'conditionParametersAV', 'conditionLabels'}));
    for structName = expStructures'
        if ~iscell(unfiltered.(structName{1})); continue; end
        filtered.(structName{1}) = filtered.(structName{1})(criterion,:);
        if length(uniquecell(filtered.(structName{1}))) == 1; filtered.(structName{1}) = filtered.(structName{1}){1}; end
    end
    if ~iscell(filtered.conditionParametersAV) && isfield(filtered, 'grids') && iscell(filtered.grids)
        filtered.grids = filtered.grids{1}; 
    elseif isfield(filtered, 'grids') && ~isstruct(filtered.grids)
        filtered.grids = filtered.grids(criterion,:);
    end
    
    filtered.index = filtIndex(filtered.index, 'expBy', criterion);
   
    expNumbers = find(criterion);
    subNumbers = unique(filtered.index.expBySubject);
    filtered = prc.filtBlock(filtered, ismember(filtered.index.trialByExp, expNumbers), 'trials');
    filtered.index = reIndex(filtered.index, 'BySubject', subNumbers);
    filtered.index = reIndex(filtered.index, 'ByExp', expNumbers);      
    
    filtered.number.trialsPerExp = arrayfun(@(x) sum(filtered.index.trialByExp==x), (1:length(expNumbers))');
    filtered.number.experiments = length(filtered.index.expBySubject);
    filtered.number.subjects = length(unique(filtered.index.expBySubject));
    
    if isfield(filtered, 'ephys') && any(filtered.index.penetrationByExp==0)
        filtered = prc.filtBlock(filtered, filtered.index.penetrationByExp~=0, 'penetrations');
    end
end

if strcmpi(filterTag, 'trials')
    trialStructures = blockFields(contains(blockFields, {'trialClass', 'timings', 'timeline', 'stim', 'inactivation', 'outcome'}));
    for structName = trialStructures'
        filtered.(structName{1}) = filterStructRows(filtered.(structName{1}), criterion);
    end
    filtered.index = filtIndex(filtered.index, 'trialBy', criterion);
    filtered.number.trials = length(filtered.index.trialByExp);
end

if strcmpi(filterTag, 'penetrations')
    filtered.ephys = filterStructRows(filtered.ephys, criterion);
    filtered.ephys.penetrationDetails = filterStructRows(filtered.ephys.penetrationDetails, criterion);
    spikeFilter = ismember(filtered.index.spikeByPenetration,find(criterion));
    clusterFilter = ismember(filtered.index.clusterByPenetration,find(criterion));
    
    filtered.ephys.spike = filterStructRows(filtered.ephys.spike, spikeFilter);
    filtered.ephys.cluster = filterStructRows(filtered.ephys.cluster, clusterFilter);
    filtered.ephys.cluster.templates = filtered.ephys.cluster.templates(criterion);
    
    filtered.index = filtIndex(filtered.index, 'penetrationBy', criterion);
    filtered.index = filtIndex(filtered.index, 'spikeBy', spikeFilter);
    filtered.index = filtIndex(filtered.index, 'clusterBy', clusterFilter);

    filtered.index = reIndex(filtered.index, 'ByPen', find(criterion));
    exp2Keep = ismember(1:filtered.number.experiments, unique(filtered.index.penetrationByExp));
    filtered = prc.filtBlock(filtered, exp2Keep(:), 'experiments');
    
    filtered.number.penetrations = sum(criterion);
    filtered.number.clusters = length(filtered.ephys.cluster.amplitudes);
    filtered.number.spikes = length(filtered.ephys.spike.amplitudes);
    filtered.number.clustersPerPenetration = arrayfun(@(x) sum(filtered.index.clusterByPenetration==x), (1:sum(criterion))');
end

end

function filtered = filterStructRows(unfiltered, criterion)
filtered = unfiltered;
fieldNames = fields(unfiltered);
for fieldName = fieldNames'
    if size(unfiltered.(fieldName{1}),1)~=length(criterion); continue; end
    filtered.(fieldName{1}) = unfiltered.(fieldName{1})(criterion,:);
end
end

function filtered = filtIndex(unfiltered, tag, criterion)
filtered = unfiltered;
indexNames = fields(unfiltered);
selectedNames = indexNames(contains(indexNames, {tag}));
for fieldName = selectedNames'
    filtered.(fieldName{1}) = filtered.(fieldName{1})(criterion); 
end
end

function filtered = reIndex(unfiltered, tag, oldNumbers)
indexNames = fields(unfiltered);
filtered = unfiltered;
selectedNames = indexNames(contains(indexNames, {tag}));
for fieldName = selectedNames'
    [~, filtered.(fieldName{1})] = ismember(filtered.(fieldName{1}), oldNumbers);
end
end
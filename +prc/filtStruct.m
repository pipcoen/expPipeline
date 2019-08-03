function filtered = filtStruct(unfiltered, criterion, filterTag)
%% Function to fulter structure.
if ~exist('criterion', 'var'); filtered = []; return; end
if ~exist('filterTag', 'var'); filterTag = 'regular'; end
filtered = unfiltered;
if contains(filterTag, {'subject'; 'experiment'})
    if iscell(criterion); experiments2keep = contains(unfiltered.subject, criterion)>0;
    else, experiments2keep = ismember(cell2mat(unfiltered.subExpPenLink(:,contains({'subject'; 'experiment'}, filterTag))), criterion)>0;
    end
    penetrations2keep = cell2mat(filtered.subExpPenLink(experiments2keep,3));
end
if strcmpi(filterTag, 'penetration')
    experiments2keep = cellfun(@(x) any(criterion==x), unfiltered.subExpPenLink(:,3));
    penetrations2keep = criterion;
end
if exist('experiments2keep', 'var')
    if length(experiments2keep) > 1; filtered = prc.filtStruct(filtered, experiments2keep); end
    filtered = prc.filtStruct(filtered, ismember(filtered.trialExperimentIdx, cell2mat(filtered.subExpPenLink(:,2)))>0);
    if isfield(filtered, 'eph_clusterPenetrationIdx')
        filtered = prc.filtStruct(filtered, ismember(filtered.eph_clusterPenetrationIdx, penetrations2keep)>0);
        filtered = prc.filtStruct(filtered, ismember(filtered.eph_spikePenetrationIdx, penetrations2keep)>0);
        filtered.subExpPenLink(:,3) = cellfun(@(x) x(ismember(x, penetrations2keep)),filtered.subExpPenLink(:,3), 'uni', 0);
    end
end

if strcmpi(filterTag, 'regular')
        criterion = criterion>0;
        refLength = length(criterion);
        fieldNames = fields(filtered);
        for fieldName = fieldNames'
            currField = fieldName{1};
            fieldSize = size(filtered.(currField));
            if ~any(fieldSize == refLength); continue; end
            filtered.(currField) = filtered.(currField)(criterion,:);
        end
end
end

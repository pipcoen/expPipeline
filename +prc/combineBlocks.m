function comBlks = combineBlocks(blks)
%% Function to combine and filter block files.
[uniqueSubjects, ~, expBySubject] = unique({blks.subject}');
if length(blks) > 255; acc1 = 'uint16'; else, acc1 = 'uint8'; end
numOfTrials = arrayfun(@(x) size(x.timings.trialStartEnd(:,1),1), blks);
trialBySubject = cell2mat(arrayfun(@(x,y) x*ones(y,1, 'uint8'), expBySubject, numOfTrials, 'uni', 0));
trialByExp = cell2mat(arrayfun(@(x,y) x*ones(y,1, acc1), cumsum(expBySubject*0+1), numOfTrials, 'uni', 0));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ephysAvailable = arrayfun(@(x) isfield(x, 'ephys'), blks);
numOfSites = arrayfun(@(x) length(x.ephys), blks(ephysAvailable));
if any(numOfSites); ephysExists = 1; else, ephysExists = 0; end

if ephysExists
    if sum(numOfSites) > 255; acc2 = 'uint16'; else, acc2 = 'uint8'; end
    ephys = vertcat(blks(numOfSites>0).ephys);
    numOfSpikes = arrayfun(@(x) length(x.spike.amplitudes),ephys);
    numOfClusters = arrayfun(@(x) length(x.cluster.amplitudes),ephys);
    numOfSpikesPerCluster = cell2mat(arrayfun(@(x) accumarray(x.spike.clusterNumber,1,[length(x.cluster.depths),1]),ephys, 'uni', 0));
    clusterNumber = cell2mat(arrayfun(@(x) (1:length(x.cluster.amplitudes))',ephys, 'uni', 0));
    
    ephys = prc.catStructs(ephys);
    ephys.penetration.subjectRef = cell2mat(arrayfun(@(x,y) x*ones(y,1), expBySubject, numOfSites, 'uni', 0));
    ephys.penetration.expRef = cell2mat(arrayfun(@(x,y) x*ones(y,1), cumsum(numOfSites*0+1), numOfSites, 'uni', 0));
    ephys.penetration.numOfClusters = numOfClusters;
    ephys.penetration.numOfSpikes = numOfSpikes;
    
    penByExp = ephys.penetration.expRef;
    ephys.spike.penetrationRef = cell2mat(arrayfun(@(x,y) x*ones(y,1, acc2), cumsum(penByExp*0+1), numOfSpikes, 'uni', 0));   
    ephys.cluster.number = clusterNumber;
    ephys.cluster.numOfSpikes = numOfSpikesPerCluster;
    ephys.cluster.subjectRef = cell2mat(arrayfun(@(x,y) x*ones(y,1, 'uint8'), ephys.penetration.subjectRef, numOfClusters, 'uni', 0));
    ephys.cluster.expRef = cell2mat(arrayfun(@(x,y) x*ones(y,1, acc1), penByExp, numOfClusters, 'uni', 0));
    ephys.cluster.penetrationRef = cell2mat(arrayfun(@(x,y) x*ones(y,1, acc2), cumsum(penByExp*0+1), numOfClusters, 'uni', 0));   
    
    expDets = kil.getExpDets(uniqueSubjects(ephys.penetration.subjectRef), {blks(penByExp).expDate}', {blks(penByExp).expNum}', ephys.penetration.folder);
    expDets = prc.catStructs(expDets);
    
    fields2copy = fields(expDets); 
    for i = fields2copy'; ephys.penetration.(i{1}) = expDets.(i{1}); end
    ephys.penetration = prc.chkThenRemoveFields(ephys.penetration, {'subject'; 'expDate'; 'expNum'; 'estDepth'});
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comBlks.tot.subjects = length(uniqueSubjects);
comBlks.tot.experiments = length(blks);
comBlks.tot.trials = sum(numOfTrials);
if ephysExists
    comBlks.tot.penetrations = sum(numOfSites);
    comBlks.tot.clusters = sum(numOfClusters);
    comBlks.tot.spikes = sum(numOfSpikes);
end

comBlks.exp.subjectRef = expBySubject;
perExpFields = {'subject', 'expDate', 'expNum', 'rigName', 'expType', 'expDef', 'conditionParametersAV', 'conditionLabels'};
if ~any(contains({blks.expDef}', 'Passive')); perExpFields = [perExpFields, 'performanceAVM']; end
for i = perExpFields; comBlks.exp.(i{1}) = {blks.(i{1})}'; end
blks = prc.chkThenRemoveFields(blks, [perExpFields, 'params', 'ephys', 'grids']);
comBlks.exp.numOfTrials = numOfTrials;

comBlks.tri.subjectRef = trialBySubject;
comBlks.tri.expRef = trialByExp;

if isfield(blks, 'timeline') && any(arrayfun(@(x) ~isempty(x.timeline), blks))
    timelineAvailable = arrayfun(@(x) isfield(x, 'timeline') & ~isempty(x.timeline), blks);
    nanTimeline = blks(find(timelineAvailable,1)).timeline;
    timelineFields = fields(nanTimeline);
    for i = timelineFields'
        if iscell(nanTimeline.(i{1})); nanTimeline.(i{1}) = {NaN};
        else, nanTimeline.(i{1}) = nanTimeline.(i{1})(1,:)*NaN;
        end
    end
    for i = find(~timelineAvailable)'
        for j = timelineFields'; blks(i).timeline.(j{1}) = repmat(nanTimeline.(j{1}), numOfTrials(i), 1); end
    end
    for i = 1:length(blks); blks(i).timeline = prc.chkThenRemoveFields(blks(i).timeline, {'alignment';'frameTimes'}); end
else, blks = prc.chkThenRemoveFields(blks, {'timeline'});
end

catBlks = prc.catStructs(blks);
trialFields = {'trialClass', 'timings', 'timeline', 'stim', 'inactivation', 'outcome'};
for i = trialFields
    if isfield(blks, i{1}); comBlks.tri.(i{1}) = catBlks.(i{1}); end
end
blks = prc.chkThenRemoveFields(blks, trialFields);

if ephysExists
    comBlks.pen = ephys.penetration;
    comBlks.clu = ephys.cluster; 
    comBlks.spk = ephys.spike; 
end

if ~isempty(fields(blks)); warning('Unexpected fields in blocks!'); end


function cellRaster(blk,eventTimes,trialGroups,sortRule)
% AP_cellraster(eventTimes,trialGroups,sortRule)
%
% Raster viewer for Neuropixels
%
% eventTimes - vector or cell array of times
% trialGroups - categorical vectors of trial types (optional)
% sortRule - sorting for scrolling through units (default = depth)
%
% These variables are required in the base workspace: 
% ephClusterTemplates (output from kilosort)
% channelPositions (output from kilosort)
% clusterDepths (calculated from center-of-mass or otherwise)
% spikeTimes_timeline (spikeTimes aligned to timeline)
% ephSpikeCluster (output from kilosort)
% spikeAmps (output from kilosort)
%
% Controls: 
% up/down - switch between units (clicking on unit also selects)
% left/right - switch between alignments (if multiple)
% u - go to unit number
% s - change session nnumber

guiData = struct;
guiData.clusterSession = blk.ephClusterSession;
guiData.clusterSite = blk.ephClusterSession;
guiData.spikeSession = blk.ephSpikeSession;
guiData.spikeSite = blk.ephSpikeSession;
guiData.blk = blk;

% Initiate eventTimes
% (align times required as input)
if ~exist('eventTimes','var'); error('No align times'); end
if exist('sortRule','var'); guiData.sortRule = sortRule; else; guiData.sortRule = 'sig'; end

% (put eventTimes into cell array if it isn't already, standardize dim)
if ~iscell(eventTimes); eventTimes = {eventTimes}; end
eventTimes = reshape(cellfun(@(x) reshape(x,[],1),eventTimes,'uni',false),1,[]);
if size(eventTimes{1},2)==1
    if length(eventTimes{1}) == length(blk.sessionIdx); eventTimes{1} = [eventTimes{1} blk.sessionIdx];
    else, error('Need to know which sessions each event relates to');
    end
end

% Initiate trialGroups
% (if no align groups specified, create one group for each alignment)
if ~exist('trialGroups','var') || isempty(trialGroups); trialGroups =  cellfun(@(x) ones(size(x,1),1),eventTimes,'uni',false); end
if ~iscell(trialGroups); trialGroups = {trialGroups}; end
trialGroups = reshape(trialGroups,1,[]);

% (replicate for each eventTimes if only one)
if length(eventTimes) > 1 && length(trialGroups) == 1
   trialGroups = repmat(trialGroups,size(eventTimes));
elseif length(eventTimes) > 1 && length(trialGroups) < length(eventTimes)
    error('Mismatching align time/group sets')
end

% (check group against time dimensions, orient align times x groups)
group_dim = cellfun(@(align,group) find(ismember(size(group),length(align))),eventTimes,trialGroups,'uni',false);
if any(cellfun(@isempty,group_dim))
    error('Mismatching times/groups within align set')
end
trialGroups = cellfun(@(groups,dim) shiftdim(groups,dim-1),trialGroups,group_dim,'uni',false);

trialGroups = {trialGroups{1}(~isnan(eventTimes{1}(:,1)),:)};
eventTimes = {eventTimes{1}(~isnan(eventTimes{1}(:,1)),:)};

% (if there isn't an all ones category first, make one)
guiData.trialGroups = cellfun(@(x) padarray(x,[0,1-all(x(:,1) == 1)],1,'pre'),trialGroups,'uni',false);

% Package gui data
cellrasterGui = figure('color','w');
guiData.currSession = 1;
guiData.currSite = 1;
guiData.eventTimes = eventTimes;

guiData.unitAxes = subplot(5,5,1:5:20,'YDir','reverse'); hold on;
xlim([-0.1,1]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

guiData.waveformAxes = subplot(5,5,2:5:20,'visible','off'); hold on;
linkaxes([guiData.unitAxes,guiData.waveformAxes],'y');

guiData.psthAxes = subplot(5,5,[3,4,5],'YAxisLocation','right'); hold on;
xlabel('Time from event (s)');
ylabel('Spikes/s/trial');

guiData.rasterAxes = subplot(5,5,[8,9,10,13,14,15,18,19,20],'YDir','reverse','YAxisLocation','right'); hold on;
xlabel('Time from event (s)');
ylabel('Trial');

guiData.amplitudeAxes = subplot(5,5,21:25); hold on;
xlabel('Experiment time (s)');
ylabel('Template amplitude');
axis tight

guidata(cellrasterGui, guiData);
cycleSession(cellrasterGui);
end

function cycleSession(cellrasterGui)
guiData = guidata(cellrasterGui);
clustIdx = guiData.blk.ephClusterSession==guiData.currSession & guiData.blk.ephClusterSite==guiData.currSite;
spikeIdx = guiData.blk.ephSpikeSession==guiData.currSession & guiData.blk.ephSpikeSite==guiData.currSite;
blk = guiData.blk;

if isfield(guiData, 'title');delete(guiData.title); end
[~, guiData.title] = plt.suplabel(sprintf('%s: %s Site %d', ...
    blk.subject{min([guiData.currSession length(blk.subject)])},blk.params(guiData.currSession).expDate,guiData.currSite),'t');
ephClusterTemplates = blk.ephClusterTemplates{guiData.currSession}{guiData.currSite};
channelPositions = blk.ephChannelMap{guiData.currSession}{guiData.currSite};
spikeTimes = blk.ephSpikeTimes(spikeIdx);
spikeCluster = blk.ephSpikeCluster(spikeIdx);
spikeAmps = blk.ephSpikeAmps(spikeIdx);
clusterDepths = blk.ephClusterDepths(clustIdx);
spikeCluster = spikeCluster-min(spikeCluster)+1;

guiData.currEventTimes = {guiData.eventTimes{1}(guiData.eventTimes{1}(:,2)==guiData.currSession,1)};
guiData.currTrialGroups = {guiData.trialGroups{1}(guiData.eventTimes{1}(:,2)==guiData.currSession,:)};

switch guiData.sortRule(1:3)
    case 'dep'
        [~, guiData.currSortRule] = sort(clusterDepths);
    case 'sig'
        blk = guiData.blk; 
        blk = prc.filtStruct(blk, clustIdx); 
        blk = prc.filtStruct(blk, spikeIdx);
        blk.ephSpikeCluster = blk.ephSpikeCluster-min(blk.ephSpikeCluster)+1;
        guiData.clusterSigLevel = kil.findResponsiveCells(blk,guiData.currEventTimes{1}, [-0.5 -0.1 0.05 0.1]);
        [~, guiData.currSortRule] = sort(guiData.clusterSigLevel);
end


% Initialize figure and axes
axes(guiData.unitAxes); cla;
ylim([-50, max(channelPositions(:,2))+50]);

% (plot unit depths by depth and relative number of spikes)
normSpikeNum = mat2gray(log(accumarray(spikeCluster,1)+1));
if length(normSpikeNum)<length(clusterDepths); normSpikeNum(end+1:length(clusterDepths)) = min(normSpikeNum); end
unitDots = plot(normSpikeNum,clusterDepths,'.k','MarkerSize',20,'ButtonDownFcn',@unitClick);
currUnitDots = plot(0,0,'.r','MarkerSize',20);

% (plot of waveform across the probe)
axes(guiData.waveformAxes); cla;
ylim([-50, max(channelPositions(:,2))+50]);
waveformLines = arrayfun(@(x) plot(guiData.waveformAxes,0,0,'k','linewidth',1),1:size(ephClusterTemplates,3));

% (smoothed psth)
axes(guiData.psthAxes); cla;
maxNumGroups = max(cell2mat(cellfun(@(x) 1+sum(diff(sort(x,1),[],1) ~= 0),guiData.currTrialGroups,'uni',false)));
psthLines = arrayfun(@(x) plot(NaN,NaN,'linewidth',2,'color','k'),1:maxNumGroups);

% (raster)
axes(guiData.rasterAxes); cla;
rasterDots = scatter(NaN,NaN,5,'k','filled');

% (spike amplitude across the recording)
axes(guiData.amplitudeAxes); cla;
amplitudePlot = plot(NaN,NaN,'.k');
amplitudeLines = arrayfun(@(x) line([0,0],ylim,'linewidth',2),1:2);

% Set default raster times
rasterWindow = [-0.5,1];
psthBinSize = 0.001;
tBins = rasterWindow(1):psthBinSize:rasterWindow(2);
t = tBins(1:end-1) + diff(tBins)./2;
useAlign = reshape(guiData.currEventTimes{1},[],1);
tPeriEvent = useAlign + tBins;
% (handle NaNs by setting rows with NaN times to 0)
tPeriEvent(any(isnan(tPeriEvent),2),:) = 0;

% Set functions for key presses
set(cellrasterGui,'KeyPressFcn',@keyPress);

% (plots)
guiData.unitDots = unitDots;
guiData.currUnitDots = currUnitDots;
guiData.waveformLines = waveformLines;
guiData.psthLines = psthLines;
guiData.rasterDots = rasterDots;
guiData.amplitudePlot = amplitudePlot;
guiData.amplitudeLines = amplitudeLines; 

% (raster times)
guiData.t = t;
guiData.tBins = tBins;
guiData.tPeriEvent = tPeriEvent;

% (spike data)
guiData.ephClusterTemplates = ephClusterTemplates;
guiData.channelPositions = channelPositions;
guiData.spikeTimes = spikeTimes;
guiData.spikeCluster = spikeCluster;
guiData.spikeAmps = spikeAmps;

% (current settings)
guiData.currUnit = guiData.currSortRule(1);
guiData.currAlign = 1;
guiData.currGroup = 1;

% Upload gui data and draw
guidata(cellrasterGui, guiData);
updatePlot(cellrasterGui);
end


function updatePlot(cellrasterGui)
% Get guidata
guiData = guidata(cellrasterGui);

% Turn on/off the appropriate graphics
set(guiData.rasterDots,'visible','on');

% Plot depth location on probe
unitX = get(guiData.unitDots,'XData');
unitY = get(guiData.unitDots,'YData');
set(guiData.currUnitDots,'XData',unitX(guiData.currUnit), 'YData',unitY(guiData.currUnit));

% Plot waveform across probe (reversed YDir, flip Y axis and plot depth)
templateXScale = 7;
templateYScale = 250;

templateY = permute(mean(guiData.ephClusterTemplates(guiData.currUnit,:,:),1),[3,2,1]);
templateY = -templateY*templateYScale + guiData.channelPositions(:,2);
templateX = (1:size(guiData.ephClusterTemplates,2)) + guiData.channelPositions(:,1)*templateXScale;

templateChannelAmp = range(guiData.ephClusterTemplates(guiData.currUnit,:,:),2);
templateThresh = max(templateChannelAmp,[],3)*0.2;
templateUseChannels = any(templateChannelAmp > templateThresh,1);
[~,maxChannel] = max(max(abs(guiData.ephClusterTemplates(guiData.currUnit,:,:)),[],2),[],3);

arrayfun(@(ch) set(guiData.waveformLines(ch),'XData',templateX(ch,:),'YData',templateY(ch,:)),1:size(guiData.ephClusterTemplates,3));
arrayfun(@(ch) set(guiData.waveformLines(ch),'Color','r'),find(templateUseChannels));
arrayfun(@(ch) set(guiData.waveformLines(ch),'Color','k'),find(~templateUseChannels));
set(guiData.waveformLines(maxChannel),'Color','b');

% Bin spikes (use only spikes within time range, big speed-up)
currSpikeIdx = ismember(guiData.spikeCluster,guiData.currUnit);
currRasterSpikeTimes = guiData.spikeTimes(currSpikeIdx);
currRasterSpikeTimes(currRasterSpikeTimes < min(guiData.tPeriEvent(:)) | ...
    currRasterSpikeTimes > max(guiData.tPeriEvent(:))) = [];

curr_raster = cell2mat(arrayfun(@(x) ...
    histcounts(currRasterSpikeTimes,guiData.tPeriEvent(x,:)), ...
    [1:size(guiData.tPeriEvent,1)]','uni',false));

% Set color scheme
currGroup = guiData.currTrialGroups{guiData.currAlign}(:,guiData.currGroup);
if length(unique(currGroup)) == 1
    % Black if one group
    group_colors = [0,0,0];
elseif length(unique(sign(currGroup(currGroup ~= 0)))) == 1
    % Black-to-red single-signed groups
    n_groups = length(unique(currGroup));
    group_colors = [linspace(0,0.8,n_groups)',zeros(n_groups,1),zeros(n_groups,1)];
elseif length(unique(sign(currGroup(currGroup ~= 0)))) == 2
    % Symmetrical blue-black-red if negative and positive groups
    n_groups_pos = length(unique(currGroup(currGroup > 0)));
    group_colors_pos = [linspace(0.3,1,n_groups_pos)',zeros(n_groups_pos,1),zeros(n_groups_pos,1)];
    
    n_groups_neg = length(unique(currGroup(currGroup < 0)));
    group_colors_neg = [zeros(n_groups_neg,1),zeros(n_groups_neg,1),linspace(0.3,1,n_groups_neg)'];
    
    n_groups_zero = length(unique(currGroup(currGroup == 0)));
    group_colors_zero = [zeros(n_groups_zero,1),zeros(n_groups_zero,1),zeros(n_groups_zero,1)];
    
    group_colors = [flipud(group_colors_neg);group_colors_zero;group_colors_pos];    
end

% Plot smoothed PSTH
smooth_size = 51;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
bin_t = mean(diff(guiData.t));

curr_psth = grpstats(curr_raster,currGroup,@(x) mean(x,1));
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0,floor(length(smWin)/2)],'replicate','both'), ...
    smWin,'valid')./bin_t;

% (set the first n lines by group, set all others to NaN)
arrayfun(@(align_group) set(guiData.psthLines(align_group), ...
    'XData',guiData.t,'YData',curr_smoothed_psth(align_group,:), ...
    'Color',group_colors(align_group,:)),1:size(curr_psth,1));
arrayfun(@(align_group) set(guiData.psthLines(align_group), ...
    'XData',NaN,'YData',NaN), ...
    size(curr_psth,1)+1:length(guiData.psthLines));

ylim(get(guiData.psthLines(1),'Parent'),[min(curr_smoothed_psth(:)), ...
    max(max(curr_smoothed_psth(:),min(curr_smoothed_psth(:))+1))]);
title(get(guiData.psthLines(1),'Parent'), ...
    ['Unit ' num2str(guiData.currUnit) ...
    ', Align ' num2str(guiData.currAlign) ...
    ', Group ' num2str(guiData.currGroup)],'FontSize',14);


% Plot raster
% (single unit mode)
[raster_y,raster_x] = find(curr_raster);
set(guiData.rasterDots,'XData',guiData.t(raster_x),'YData',raster_y);
xlim(get(guiData.rasterDots,'Parent'),[guiData.tBins(1),guiData.tBins(end)]);
ylim(get(guiData.rasterDots,'Parent'),[0,size(guiData.tPeriEvent,1)]);
% (set dot color by group)
[~,~,row_group] = unique(guiData.currTrialGroups{guiData.currAlign}(:,guiData.currGroup),'sorted');
psth_colors = get(guiData.psthLines,'color');
if iscell(psth_colors); psth_colors = cell2mat(psth_colors); end
raster_dot_color = psth_colors(row_group(raster_y),:);
set(guiData.rasterDots,'CData',raster_dot_color);


% Plot template amplitude over whole experiment
set(guiData.amplitudePlot,'XData', ...
    guiData.spikeTimes(currSpikeIdx), ...
    'YData',guiData.spikeAmps(currSpikeIdx),'linestyle','none');

[ymin,ymax] = bounds(get(guiData.amplitudePlot,'YData'));
set(guiData.amplitudeLines(1),'XData',repmat(min(guiData.tPeriEvent(:)),2,1),'YData',[ymin,ymax]);
set(guiData.amplitudeLines(2),'XData',repmat(max(guiData.tPeriEvent(:)),2,1),'YData',[ymin,ymax]);
assignin('base','guiData',guiData)
end


function keyPress(cellrasterGui,eventdata)

% Get guidata
guiData = guidata(cellrasterGui);

switch eventdata.Key
    case 's'
        numSessionSites = length(guiData.blk.ephClusterTemplates{guiData.currSession});
        if guiData.currSite == numSessionSites; guiData.currSite = 1; guiData.currSession = guiData.currSession+1;
        else; guiData.currSite = guiData.currSite+1; 
        end
        if guiData.currSession>max(guiData.clusterSession); guiData.currSession = 1; end
        guidata(cellrasterGui,guiData);
        cycleSession(cellrasterGui);
        guiData = guidata(cellrasterGui);
    
    case 'downarrow'
        % Next unit
        currUnitIdx = guiData.currUnit(1) == guiData.currSortRule;
        newUnit = guiData.currSortRule(circshift(currUnitIdx,1));
        guiData.currUnit = newUnit;
        
    case 'uparrow'
        % Previous unit
        currUnitIdx = guiData.currUnit(end) == guiData.currSortRule;
        newUnit = guiData.currSortRule(circshift(currUnitIdx,-1));
        guiData.currUnit = newUnit;
        
    case 'rightarrow'
        % Next alignment
        newAlign = guiData.currAlign + 1;
        if newAlign > length(guiData.currEventTimes)
            newAlign = 1;
        end
        useAlign = reshape(guiData.currEventTimes{newAlign},[],1);
        tPeriEvent = useAlign + guiData.tBins;
        
        % (handle NaNs by setting rows with NaN times to 0)
        tPeriEvent(any(isnan(tPeriEvent),2),:) = 0;
        
        guiData.currAlign = newAlign;
        guiData.tPeriEvent = tPeriEvent;
        guiData.currGroup = 1;
        
    case 'leftarrow'
        % Previous alignment
        newAlign = guiData.currAlign - 1;
        if newAlign < 1
            newAlign = length(guiData.currEventTimes);
        end
        useAlign = reshape(guiData.currEventTimes{newAlign},[],1);
        tPeriEvent = useAlign + guiData.tBins;
        
        % (handle NaNs by setting rows with NaN times to 0)
        tPeriEvent(any(isnan(tPeriEvent),2),:) = 0;
        
        guiData.currAlign = newAlign;
        guiData.tPeriEvent = tPeriEvent;
        guiData.currGroup = 1;

    case 'pagedown'
        % Next group
        nextGroup = guiData.currGroup + 1;
        if nextGroup > size(guiData.currTrialGroups{guiData.currAlign},2)
            nextGroup = 1;
        end
        guiData.currGroup = nextGroup;
        
    case 'pageup'
        % Previous group
        nextGroup = guiData.currGroup - 1;
        if nextGroup < 1
            nextGroup = size(guiData.currTrialGroups{guiData.currAlign},2);
        end
        guiData.currGroup = nextGroup;     
        
    case 'u'
        % Enter and go to unit
        newUnit = str2num(cell2mat(inputdlg('Go to unit:')));
        if ~ismember(newUnit,unique(guiData.spikeCluster))
            error(['Unit ' num2str(newUnit) ' not present'])
        end
        guiData.currUnit = newUnit;
        
end

% Upload gui data and draw
guidata(cellrasterGui,guiData);
updatePlot(cellrasterGui);
end

function unitClick(cellrasterGui,eventdata)

% Get guidata
guiData = guidata(cellrasterGui);

% Get the clicked unit, update current unit
unitX = get(guiData.unitDots,'XData');
unitY = get(guiData.unitDots,'YData');

[~,clicked_unit] = min(sqrt(sum(([unitX;unitY] - ...
    eventdata.IntersectionPoint(1:2)').^2,1)));

guiData.currUnit = clicked_unit;

% Upload gui data and draw
guidata(cellrasterGui,guiData);
updatePlot(cellrasterGui);

end





















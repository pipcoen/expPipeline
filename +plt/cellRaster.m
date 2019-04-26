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
% ephWaveforms (output from kilosort)
% channelPositions (output from kilosort)
% clusterDepths (calculated from center-of-mass or otherwise)
% spikeTimes_timeline (spikeTimes aligned to timeline)
% ephClusterID (output from kilosort)
% spikeAmps (output from kilosort)
%
% Controls: 
% up/down - switch between units (clicking on unit also selects)
% left/right - switch between alignments (if multiple)
% m - select depth range to plot multiunit
% u - go to unit number
% s - change session nnumber

guiData = struct;
if ~isfield(blk, 'ephClustSessionIdx'); guiData.clustSessions = blk.ephClusterAmps*0+1; else, guiData.clustSessions = blk.ephClustSessionIdx; end
if ~isfield(blk, 'ephSessionIdx'); guiData.spikeSessions = blk.ephSessionIdx*0+1; else, guiData.spikeSessions = blk.ephSessionIdx; end
guiData.blk = blk;

% Initiate eventTimes
% (align times required as input)
if ~exist('eventTimes','var'); error('No align times'); end
if size(eventTimes,2)==1
    if length(eventTimes) == lenghth(blk.sessionIdx); guiData.eventTimes = [guiData.eventTimes blk.sessionIdx];
    else, error('Need to know which sessions each event relates to');
    end
end
if exist('sortRule','var'); guiData.sortRule = sortRule; end

% (put eventTimes into cell array if it isn't already, standardize dim)
if ~iscell(eventTimes); eventTimes = {eventTimes}; end
eventTimes = reshape(cellfun(@(x) reshape(x,[],1),eventTimes,'uni',false),1,[]);

% Initiate trialGroups
% (if no align groups specified, create one group for each alignment)
if ~exist('trialGroups','var') || isempty(trialGroups)
   guiData.trialGroups =  cellfun(@(x) ones(size(x,1),1),eventTimes,'uni',false);
else guiData.trialGroups = trialGroups;
end

% (put groups into cell array if it isn't already)
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

% (if there isn't an all ones category first, make one)
trialGroups = cellfun(@(x) padarray(x,[0,1-all(x(:,1) == 1)],1,'pre'),trialGroups,'uni',false);

% Package gui data
cellrasterGui = figure('color','w');
guiData.selectedSession = 0;
guiData.evenTimes = eventTimes;
guidata(cellrasterGui, guiData);
cycleSession(cellrasterGui);
end

function cycleSession(cellrasterGui)
guiData = guidata(cellrasterGui);
guiData.selectedSession = mod(guiData.selectedSession+1, max(guiData.clustSessions));
% Pull standard ephys variables from base workspace
ephWaveforms = guiData.blk.ephWaveforms(guiData.clustSessions==guiData.selectedSession);
channelPositions = guiData.blk.ephChannelMap{guiData.selectedSession};
spikeTimes = guiData.blk.ephSpikeTimes(guiData.spikeSessions==guiData.selectedSession);
ephClusterID = guiData.blk.ephClusterID(guiData.spikeSessions==guiData.selectedSession);
spikeAmps = guiData.blk.ephSpikeAmps(guiData.clustSessions==guiData.selectedSession);
clusterDepths = guiData.blk.ephClusterDepths(guiData.clustSessions==guiData.selectedSession);
% clusterDepths = guiData.blk.ephClusterDepths(unique(ephClusterID));

% Sort the units by depth if not specified
if ~isfield(guiData, 'sortRule'); [~,sortRule] = sort(clusterDepths); end

% Initialize figure and axes

% (plot unit depths by depth and relative number of spikes)
unit_axes = subplot(5,5,[1:5:20],'YDir','reverse');
hold on;

norm_spike_n = mat2gray(log(accumarray(ephClusterID,1)+1));
unit_dots = plot(norm_spike_n,clusterDepths,'.k','MarkerSize',20,'ButtonDownFcn',@unit_click);
curr_unit_dots = plot(0,0,'.r','MarkerSize',20);
multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
xlim(unit_axes,[-0.1,1]);
ylim([-50, max(channelPositions(:,2))+50]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

% (plot of waveform across the probe)
waveform_axes = subplot(5,5,[2:5:20],'visible','off','YDir','reverse');
hold on;
ylim([-50, max(channelPositions(:,2))+50]);
waveform_lines = arrayfun(@(x) plot(waveform_axes,0,0,'k','linewidth',1),1:size(ephWaveforms,3));

linkaxes([unit_axes,waveform_axes],'y');

% (smoothed psth)
psth_axes = subplot(5,5,[3,4,5],'YAxisLocation','right');
hold on;
max_n_groups = max(cell2mat(cellfun(@(x) 1+sum(diff(sort(x,1),[],1) ~= 0),trialGroups,'uni',false)));
psth_lines = arrayfun(@(x) plot(NaN,NaN,'linewidth',2,'color','k'),1:max_n_groups);
xlabel('Time from event (s)');
ylabel('Spikes/s/trial');

% (raster)
raster_axes = subplot(5,5,[8,9,10,13,14,15,18,19,20],'YDir','reverse','YAxisLocation','right');
hold on;
raster_dots = scatter(NaN,NaN,5,'k','filled');
raster_image = imagesc(NaN,'visible','off'); colormap(raster_axes,hot);
xlabel('Time from event (s)');
ylabel('Trial');

% (spike amplitude across the recording)
amplitude_axes = subplot(5,5,21:25); hold on;
amplitude_plot = plot(NaN,NaN,'.k');
amplitude_lines = arrayfun(@(x) line([0,0],ylim,'linewidth',2),1:2);
xlabel('Experiment time (s)');
ylabel('Template amplitude');
axis tight

% Set default raster times
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = t_bins(1:end-1) + diff(t_bins)./2;
use_align = reshape(eventTimes{1},[],1);
t_peri_event = use_align + t_bins;
% (handle NaNs by setting rows with NaN times to 0)
t_peri_event(any(isnan(t_peri_event),2),:) = 0;

% Set functions for key presses
set(cellrasterGui,'KeyPressFcn',@key_press);

% (plots)
guiData.unit_dots = unit_dots;
guiData.curr_unit_dots = curr_unit_dots;
guiData.multiunit_lines = multiunit_lines;
guiData.waveform_lines = waveform_lines;
guiData.psth_lines = psth_lines;
guiData.raster_dots = raster_dots;
guiData.raster_image = raster_image;
guiData.amplitude_plot = amplitude_plot;
guiData.amplitude_lines = amplitude_lines; 

% (raster times)
guiData.t = t;
guiData.t_bins = t_bins;
guiData.t_peri_event = t_peri_event;

% (user inputs)
guiData.eventTimes = eventTimes;
guiData.trialGroups = trialGroups;
guiData.sortRule = sortRule;

% (spike data)
guiData.ephWaveforms = ephWaveforms;
guiData.channelPositions = channelPositions;
guiData.spikeTimes = spikeTimes;
guiData.ephClusterID = ephClusterID;
guiData.spikeAmps = spikeAmps;

% (current settings)
guiData.curr_unit = sortRule(1);
guiData.curr_align = 1;
guiData.curr_group = 1;

% Upload gui data and draw
guidata(cellrasterGui, guiData);
update_plot(cellrasterGui);
end


function update_plot(cellrasterGui,eventdata)

% Get guidata
guiData = guidata(cellrasterGui);

% Turn on/off the appropriate graphics
if length(guiData.curr_unit) == 1
    set(guiData.raster_dots,'visible','on');
    set(guiData.multiunit_lines,'visible','off');
    set(guiData.raster_image,'visible','off');
elseif length(guiData.curr_unit) > 1
    set(guiData.raster_dots,'visible','off');
    set(guiData.multiunit_lines,'visible','on');
    set(guiData.raster_image,'visible','on');
end

% Plot depth location on probe
unit_x = get(guiData.unit_dots,'XData');
unit_y = get(guiData.unit_dots,'YData');
set(guiData.curr_unit_dots,'XData',unit_x(guiData.curr_unit), ...
    'YData',unit_y(guiData.curr_unit));

% Plot waveform across probe (reversed YDir, flip Y axis and plot depth)
template_xscale = 7;
template_yscale = 5;

template_y = permute(mean(guiData.ephWaveforms(guiData.curr_unit,:,:),1),[3,2,1]);
template_y = -template_y*template_yscale + guiData.channelPositions(:,2);
template_x = (1:size(guiData.ephWaveforms,2)) + guiData.channelPositions(:,1)*template_xscale;

template_channel_amp = range(guiData.ephWaveforms(guiData.curr_unit,:,:),2);
template_thresh = max(template_channel_amp,[],3)*0.2;
template_use_channels = any(template_channel_amp > template_thresh,1);
[~,max_channel] = max(max(abs(guiData.ephWaveforms(guiData.curr_unit,:,:)),[],2),[],3);

arrayfun(@(ch) set(guiData.waveform_lines(ch),'XData',template_x(ch,:),'YData',template_y(ch,:)),1:size(guiData.ephWaveforms,3));
arrayfun(@(ch) set(guiData.waveform_lines(ch),'Color','r'),find(template_use_channels));
arrayfun(@(ch) set(guiData.waveform_lines(ch),'Color','k'),find(~template_use_channels));
set(guiData.waveform_lines(max_channel),'Color','b');

% Bin spikes (use only spikes within time range, big speed-up)
curr_spikes_idx = ismember(guiData.ephClusterID,guiData.curr_unit);
curr_raster_spikeTimes = guiData.spikeTimes(curr_spikes_idx);
curr_raster_spikeTimes(curr_raster_spikeTimes < min(guiData.t_peri_event(:)) | ...
    curr_raster_spikeTimes > max(guiData.t_peri_event(:))) = [];

curr_raster = cell2mat(arrayfun(@(x) ...
    histcounts(curr_raster_spikeTimes,guiData.t_peri_event(x,:)), ...
    [1:size(guiData.t_peri_event,1)]','uni',false));

% Set color scheme
curr_group = guiData.trialGroups{guiData.curr_align}(:,guiData.curr_group);
if length(unique(curr_group)) == 1
    % Black if one group
    group_colors = [0,0,0];
elseif length(unique(sign(curr_group(curr_group ~= 0)))) == 1
    % Black-to-red single-signed groups
    n_groups = length(unique(curr_group));
    group_colors = [linspace(0,0.8,n_groups)',zeros(n_groups,1),zeros(n_groups,1)];
elseif length(unique(sign(curr_group(curr_group ~= 0)))) == 2
    % Symmetrical blue-black-red if negative and positive groups
    n_groups_pos = length(unique(curr_group(curr_group > 0)));
    group_colors_pos = [linspace(0.3,1,n_groups_pos)',zeros(n_groups_pos,1),zeros(n_groups_pos,1)];
    
    n_groups_neg = length(unique(curr_group(curr_group < 0)));
    group_colors_neg = [zeros(n_groups_neg,1),zeros(n_groups_neg,1),linspace(0.3,1,n_groups_neg)'];
    
    n_groups_zero = length(unique(curr_group(curr_group == 0)));
    group_colors_zero = [zeros(n_groups_zero,1),zeros(n_groups_zero,1),zeros(n_groups_zero,1)];
    
    group_colors = [flipud(group_colors_neg);group_colors_zero;group_colors_pos];    
end

% Plot smoothed PSTH
smooth_size = 51;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
bin_t = mean(diff(guiData.t));

curr_psth = grpstats(curr_raster,curr_group,@(x) mean(x,1));
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0,floor(length(smWin)/2)],'replicate','both'), ...
    smWin,'valid')./bin_t;

% (set the first n lines by group, set all others to NaN)
arrayfun(@(align_group) set(guiData.psth_lines(align_group), ...
    'XData',guiData.t,'YData',curr_smoothed_psth(align_group,:), ...
    'Color',group_colors(align_group,:)),1:size(curr_psth,1));
arrayfun(@(align_group) set(guiData.psth_lines(align_group), ...
    'XData',NaN,'YData',NaN), ...
    size(curr_psth,1)+1:length(guiData.psth_lines));

ylim(get(guiData.psth_lines(1),'Parent'),[min(curr_smoothed_psth(:)), ...
    max(max(curr_smoothed_psth(:),min(curr_smoothed_psth(:))+1))]);
if length(guiData.curr_unit) == 1
title(get(guiData.psth_lines(1),'Parent'), ...
    ['Unit ' num2str(guiData.curr_unit) ...
    ', Align ' num2str(guiData.curr_align) ...
    ', Group ' num2str(guiData.curr_group)],'FontSize',14);
elseif length(guiData.curr_unit) > 1
title(get(guiData.psth_lines(1),'Parent'), ...
    ['Multiunit, Align ' num2str(guiData.curr_align)],'FontSize',14);
end

% Plot raster
if length(guiData.curr_unit) == 1
    % (single unit mode)
    [raster_y,raster_x] = find(curr_raster);
    set(guiData.raster_dots,'XData',guiData.t(raster_x),'YData',raster_y);
    xlim(get(guiData.raster_dots,'Parent'),[guiData.t_bins(1),guiData.t_bins(end)]);
    ylim(get(guiData.raster_dots,'Parent'),[0,size(guiData.t_peri_event,1)]);
    % (set dot color by group)
    [~,~,row_group] = unique(guiData.trialGroups{guiData.curr_align}(:,guiData.curr_group),'sorted');
    psth_colors = get(guiData.psth_lines,'color');
    if iscell(psth_colors); psth_colors = cell2mat(psth_colors); end
    raster_dot_color = psth_colors(row_group(raster_y),:);
    set(guiData.raster_dots,'CData',raster_dot_color);
    
elseif length(guiData.curr_unit) > 1
    % (multiunit mode)
    raster_heatmap = imgaussfilt(curr_raster,[5,10]);
    set(guiData.raster_image,'XData',guiData.t,'YData', ...
        1:size(guiData.t_peri_event,1),'CData',raster_heatmap);
    caxis(get(guiData.raster_image,'Parent'),prctile(raster_heatmap(:),[0.05,99.5]));
    
end

% Plot template amplitude over whole experiment
if length(guiData.curr_unit) == 1
    set(guiData.amplitude_plot,'XData', ...
        guiData.spikeTimes(curr_spikes_idx), ...
        'YData',guiData.spikeAmps(curr_spikes_idx),'linestyle','none');
elseif length(guiData.curr_unit) > 1
    long_bin_size = 60;
    long_bins = guiData.spikeTimes(1):long_bin_size:guiData.spikeTimes(end);
    long_bins_t = long_bins(1:end-1) + diff(long_bins)/2;
    long_spikes_binned = discretize(guiData.spikeTimes,long_bins);
    amplitude_binned = accumarray(long_spikes_binned(curr_spikes_idx & ~isnan(long_spikes_binned)), ...
        guiData.spikeAmps(curr_spikes_idx & ~isnan(long_spikes_binned)),size(long_bins_t'),@nansum,NaN);
    set(guiData.amplitude_plot,'XData',long_bins_t,'YData',amplitude_binned,'linestyle','-');
end

[ymin,ymax] = bounds(get(guiData.amplitude_plot,'YData'));
set(guiData.amplitude_lines(1),'XData',repmat(min(guiData.t_peri_event(:)),2,1),'YData',[ymin,ymax]);
set(guiData.amplitude_lines(2),'XData',repmat(max(guiData.t_peri_event(:)),2,1),'YData',[ymin,ymax]);

end


function key_press(cellrasterGui,eventdata)

% Get guidata
guiData = guidata(cellrasterGui);

switch eventdata.Key
    case 'downarrow'
        % Next unit
        curr_unit_idx = guiData.curr_unit(1) == guiData.sortRule;
        new_unit = guiData.sortRule(circshift(curr_unit_idx,1));
        guiData.curr_unit = new_unit;
        
    case 'uparrow'
        % Previous unit
        curr_unit_idx = guiData.curr_unit(end) == guiData.sortRule;
        new_unit = guiData.sortRule(circshift(curr_unit_idx,-1));
        guiData.curr_unit = new_unit;
        
    case 'rightarrow'
        % Next alignment
        new_align = guiData.curr_align + 1;
        if new_align > length(guiData.eventTimes)
            new_align = 1;
        end
        use_align = reshape(guiData.eventTimes{new_align},[],1);
        t_peri_event = use_align + guiData.t_bins;
        
        % (handle NaNs by setting rows with NaN times to 0)
        t_peri_event(any(isnan(t_peri_event),2),:) = 0;
        
        guiData.curr_align = new_align;
        guiData.t_peri_event = t_peri_event;
        guiData.curr_group = 1;
        
    case 'leftarrow'
        % Previous alignment
        new_align = guiData.curr_align - 1;
        if new_align < 1
            new_align = length(guiData.eventTimes);
        end
        use_align = reshape(guiData.eventTimes{new_align},[],1);
        t_peri_event = use_align + guiData.t_bins;
        
        % (handle NaNs by setting rows with NaN times to 0)
        t_peri_event(any(isnan(t_peri_event),2),:) = 0;
        
        guiData.curr_align = new_align;
        guiData.t_peri_event = t_peri_event;
        guiData.curr_group = 1;

    case 'pagedown'
        % Next group
        next_group = guiData.curr_group + 1;
        if next_group > size(guiData.trialGroups{guiData.curr_align},2)
            next_group = 1;
        end
        guiData.curr_group = next_group;
        
    case 'pageup'
        % Previous group
        next_group = guiData.curr_group - 1;
        if next_group < 1
            next_group = size(guiData.trialGroups{guiData.curr_align},2);
        end
        guiData.curr_group = next_group;
                
    case 'm'
        % Multiunit (select on unit depth plot)
        [~,multiunit_top] = ginput(1);
        set(guiData.multiunit_lines(1),'visible','on','YData',repmat(multiunit_top,1,2));
        [~,multiunit_bottom] = ginput(1);
        set(guiData.multiunit_lines(2),'visible','on','YData',repmat(multiunit_bottom,1,2));
        
        clusterDepths = get(guiData.unit_dots,'YData');
        guiData.curr_unit = find(clusterDepths >= multiunit_top & ...
            clusterDepths <= multiunit_bottom);       
        
    case 'u'
        % Enter and go to unit
        new_unit = str2num(cell2mat(inputdlg('Go to unit:')));
        if ~ismember(new_unit,unique(guiData.ephClusterID))
            error(['Unit ' num2str(new_unit) ' not present'])
        end
        guiData.curr_unit = new_unit;
        
end

% Upload gui data and draw
guidata(cellrasterGui,guiData);
update_plot(cellrasterGui);
        
end

function unit_click(cellrasterGui,eventdata)

% Get guidata
guiData = guidata(cellrasterGui);

% Get the clicked unit, update current unit
unit_x = get(guiData.unit_dots,'XData');
unit_y = get(guiData.unit_dots,'YData');

[~,clicked_unit] = min(sqrt(sum(([unit_x;unit_y] - ...
    eventdata.IntersectionPoint(1:2)').^2,1)));

guiData.curr_unit = clicked_unit;

% Upload gui data and draw
guidata(cellrasterGui,guiData);
update_plot(cellrasterGui);

end





















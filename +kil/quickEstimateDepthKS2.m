function quickEstimateDepthKS2(lfpPowerSpectra, sName)
%%
guiData = struct;
depthEstimateGUI = figure('color','w');
set(depthEstimateGUI,'KeyPressFcn',@keyPress);
guiData.imAx = axes(depthEstimateGUI);
set(guiData.imAx, 'xlim', [0 200]);
set(depthEstimateGUI, 'position', get(depthEstimateGUI, 'position').*[1 0.2 1 2])

guiData.numOfShanks = length(unique(floor(lfpPowerSpectra.xCoord/100)));
guiData.lfpPowerSpectra = lfpPowerSpectra;
guiData.currShank = 0;
if ~isfield(guiData.lfpPowerSpectra, 'surfaceEst')
    guiData.lfpPowerSpectra.surfaceEst = zeros(guiData.numOfShanks,1);
end

guiData.lfpSpectra = rand(384,200);
guiData.sName = sName;
guidata(depthEstimateGUI, guiData);
guidata(depthEstimateGUI);
updatePlot(depthEstimateGUI);
end

function updatePlot(depthEstimateGUI)
guiData = guidata(depthEstimateGUI);
currIdx = floor((guiData.lfpPowerSpectra.xCoord)/200)==guiData.currShank;

LFP = guiData.lfpPowerSpectra.powerSpectra(currIdx,:);
yCoord = guiData.lfpPowerSpectra.yCoord(currIdx);

if guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1) == 0
    guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1) = median(yCoord);
end

[~,sIdx] = sort(yCoord);
axes(guiData.imAx);
cla;
LFP = (log(LFP));
% channelMean = mean(LFP,2);
% % channels2Remove = abs(channelMean-mean(channelMean))>15;
% % LFP(channels2Remove,:) = median(LFP(:));
imagesc(1:200, yCoord(sIdx), (LFP(sIdx,:))); axis xy
hold on;
plot(xlim, guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1).*[1 1], '--k', 'linewidth', 2)
guidata(depthEstimateGUI, guiData);
end

function keyPress(depthEstimateGUI,eventdata)
guiData = guidata(depthEstimateGUI);
switch eventdata.Key
    case 'downarrow'
        % Next unit
        guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1) = guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1)-15;        
    case 'uparrow'
        % Previous unit
        guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1) = guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1)+15;
    case 'pageup'
        % Next unit
        guiData.currShank = min([guiData.numOfShanks-1, guiData.currShank+1]);
    case 'pagedown'
        % Previous unit
        guiData.currShank = max([0 guiData.currShank-1]);
    case 'n'
        % Previous unit
        guiData.lfpPowerSpectra.surfaceEst(guiData.currShank+1) = nan;
    case 's'
        disp('Saving LFP Data');
        lfpPowerSpectra = guiData.lfpPowerSpectra;
        save(guiData.sName, 'lfpPowerSpectra')
    case 'q'
        close(depthEstimateGUI);
        return;
        % Next alignment
end
guidata(depthEstimateGUI, guiData);
updatePlot(depthEstimateGUI);
end
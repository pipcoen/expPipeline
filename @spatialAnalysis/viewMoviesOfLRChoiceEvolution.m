function viewMoviesOfLRChoiceEvolution(obj)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices
%	'rea'--------------------timeToFirstMove vs fration of rightward choices
figure;
axesOpt.totalNumOfAxes = length(obj.blks);
axesOpt.btlrMargins = [80 100 80 40];
axesOpt.gapBetweenAxes = [100 0];
axesOpt.numOfRows = 3;%ceil(length(obj.blks)/5);
axesOpt.axesSize = [400 450];

for i  = 1:length(obj.blks)
    bigBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    bigBlk = prc.filtBlock(bigBlk, ~isnan(bigBlk.tri.outcome.timeToFirstMove) & bigBlk.tri.outcome.timeToFirstMove<0.5);
    contast2Use = max(abs(bigBlk.tri.stim.visContrast(bigBlk.tri.trialType.visual)));
    bigBlk = prc.filtBlock(bigBlk, (abs(bigBlk.tri.stim.visContrast) == contast2Use) | bigBlk.tri.trialType.auditory);
    
    axesHand{i} = plt.getAxes(axesOpt, i);
    title(sprintf('%s: %d Tri', cell2mat(unique(bigBlk.exp.subject)'), bigBlk.tot.trials))
    box off;
    hold on;
    obj.blks(i) = bigBlk;
    axis square
    axis off
end
figureSize = get(gcf, 'position');
mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
plt.suplabel('\fontsize{20} Normalized velocity', 'y', mainAxes);
plt.suplabel('\fontsize{20} Normalized displacement', 'x', mainAxes);


evalPnts = 0.030:0.002:0.2;
tStep = evalPnts(2)-evalPnts(1);
moveDatAVC = cell(3,2);
colorBuff = floor(length(evalPnts)/20)*2+1;
for i  = 1:length(obj.blks)
    bigBlk = obj.blks(i);
    rawWheelTV = bigBlk.tri.raw.wheelTimeValue;
    
    alignTimes = repmat({evalPnts}, bigBlk.tot.trials,1);
    bigBlk.tri.wheelPos = cellfun(@(x,y) interp1(x(:,1), x(:,2), [y y(end)+tStep], 'nearest', 'extrap')*-1, rawWheelTV, alignTimes, 'uni', 0);
    bigBlk.tri.wheelVel = cellfun(@diff, bigBlk.tri.wheelPos, 'uni', 0);
    bigBlk.tri.wheelPos = cellfun(@(x) x(1:end-1), bigBlk.tri.wheelPos, 'uni', 0);
    
    threshMoveTime = bigBlk.tri.outcome.threshMoveTime;
    allWheelPosRaw = cell2mat(arrayfun(@(x,y,z) [x{1}(z{1}<=y) nan*x{1}(z{1}>y)], bigBlk.tri.wheelPos, threshMoveTime, alignTimes, 'uni', 0));
    allWheelVelRaw = cell2mat(arrayfun(@(x,y,z) [x{1}(z{1}<=y) nan*x{1}(z{1}>y)], bigBlk.tri.wheelVel, threshMoveTime, alignTimes,  'uni', 0));
    
    allWheelPos{i,1} = normalize(allWheelPosRaw, 'scale');
    allWheelVel{i,1} = normalize(allWheelVelRaw, 'scale');
    zeroIdx{i,1} = allWheelPosRaw==0 & allWheelVelRaw==0;
    idxR{i,1} = bigBlk.tri.outcome.threshMoveDirection==2;
    idxL{i,1} = bigBlk.tri.outcome.threshMoveDirection==1;
end
% %%
% vid = VideoWriter('test2');
% vid.FrameRate = 5;
% open(vid)
% axlim = 20;
% for j = 1:length(evalPnts)
%     for i = 1:length(obj.blks)
%         axes(axesHand{i});
%         cla
%         zIdx = ~(zeroIdx{i}(:,j)>100);
%         posDat = allWheelPos{i}(zIdx,j);
%         velDat = allWheelVel{i}(zIdx,j);
%         iR = idxR{i}(zIdx);
%         iL = idxL{i}(zIdx);
%         
%         scatter(posDat(iR), velDat(iR), 'r', 'filled', 'MarkerFaceAlpha', 0.15)
%         hold on
%         scatter(posDat(iL), velDat(iL), 'b', 'filled', 'MarkerFaceAlpha', 0.15)
%         xlim([-1 1]*axlim);
%         ylim([-1 1]*axlim);
%     end
%     axTitle = plt.suplabel(['\fontsize{20} Time from stimulus onset: ' num2str(evalPnts(j)) 'ms'], 't', mainAxes);
%     writeVideo(vid,getframe(gcf))
%     delete(axTitle);
% end
% close(vid)


%%
vid = VideoWriter('test2');
vid.FrameRate = 5;
open(vid)
cLevels = 10;
bWidth = 0.5;
colorMap = sqrt(plt.redBlueMap(cLevels*2+colorBuff));
colL = flipud(colorMap(1:cLevels,:));
colR = colorMap(end-cLevels+1:end,:);
cPnts = -2:0.1:2;
[xDat,yDat] = meshgrid(cPnts', cPnts');
for j = 1:length(evalPnts)
    for i = 1:length(obj.blks)
        zIdx = ~(zeroIdx{i}(:,j)>100);
        posDat = allWheelPos{i}(zIdx,j);
        velDat = allWheelVel{i}(zIdx,j);
        iR = idxR{i}(zIdx);
        iL = idxL{i}(zIdx);
        
        if all(isnan(velDat)); continue; end
        if all(isnan(posDat)); continue; end
        
        
        axes(axesHand{i});
        cla
        plot([0 0], [min(cPnts), max(cPnts)], '--k', 'linewidth', 2);
        plot([min(cPnts), max(cPnts)], [0 0], '--k', 'linewidth', 2);
%         if sum(~zeroIdx{i}(:,j)) <500; continue; end
       
        
        [fDat] = ksdensity([posDat(iR), velDat(iR)], [xDat(:), yDat(:)], 'bandwidth', bWidth);
        contMat = contourc(cPnts,cPnts, reshape(fDat, [length(cPnts),length(cPnts)]), cLevels);
        [x,y] = C2xyz(contMat);
        for k = cLevels-3:cLevels
            plot(x{k}, y{k}, 'color', colR(k,:))
        end
        
        hold on
        [fDat] = ksdensity([posDat(iL), velDat(iL)], [xDat(:), yDat(:)], 'bandwidth', bWidth);
        contMat = contourc(cPnts,cPnts, reshape(fDat, [length(cPnts),length(cPnts)]), cLevels);
        [x,y] = C2xyz(contMat);
        for k = cLevels-5:cLevels
            plot(x{k}, y{k}, 'color', colL(k,:))
        end
        xlim([min(cPnts), max(cPnts)]);
        ylim([min(cPnts), max(cPnts)]);
    end
    axTitle = plt.suplabel(['\fontsize{20} Time from stimulus onset: ' num2str(evalPnts(j)) 'ms'], 't', mainAxes);
    writeVideo(vid,getframe(gcf))
    delete(axTitle);
end
close(vid)
%%
xlim([-1 1]);
ylim([-1 1]);
axis square
plot([-1 1], [0 0], '--k', 'linewidth', 3)
plot([0 0], [-1 1], '--k', 'linewidth', 3)
xlabel('Nomalized Displacement');
ylabel('Nomalized Velocity');


end

function [x,y,z] = C2xyz(C)
m(1)=1;
n=1;
try
    while n<length(C)
        n=n+1;
        m(n) = m(n-1)+C(2,m(n-1))+1;
        
    end
end
for nn = 1:n-2
    x{nn} = C(1,m(nn)+1:m(nn+1)-1);
    y{nn} = C(2,m(nn)+1:m(nn+1)-1);
    if nargout==3
        z(nn) = C(1,m(nn));
    end
end
end


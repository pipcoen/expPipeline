function fastReactionTimes
%%
s = spatialAnalysis('all', 'm2ephysgood', 1, 1);
%%
blk = s.blks;
timelinePaths = cell(blk.tot.experiments,1);
for j = 1:blk.tot.experiments
    pathInfo.subject = blk.exp.subject{j};
    pathInfo.expDate = blk.exp.expDate{j};
    pathInfo.expNum = blk.exp.expNum{j};
    timelinePaths{j} = prc.pathFinder('servertimeline', pathInfo);
end
firstMoveTimes = blk.tri.timeline.firstMoveTimes;
timelineRefs = timelinePaths(blk.tri.expRef);
reactionTimes = firstMoveTimes - nanmin([blk.tri.timeline.audStimPeriodOnOff(:,1) blk.tri.timeline.visStimPeriodOnOff(:,1)], [], 2);

times2LookAt = find(reactionTimes>0 & reactionTimes<0.01);
sessions2LookAt = unique(blk.tri.expRef(times2LookAt));

%%

currSession = 3;%mode(blk.tri.expRef(times2LookAt));
sessionTimes = firstMoveTimes(times2LookAt(blk.tri.expRef(times2LookAt)==currSession));
load(timelinePaths{currSession});
timeline = Timeline;
inputNames = {timeline.hw.inputs.name}';                      %List of inputs to the timeline file
timelineTime = timeline.rawDAQTimestamps;                     %Timestamps in the timeline file
timelinehWeelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
timelinehWeelPosition(timelinehWeelPosition > 2^31) = timelinehWeelPosition(timelinehWeelPosition > 2^31) - 2^32;
cla;
plot(timelineTime, timelinehWeelPosition)
hold on;

sessionTimeIdx = round(sessionTimes*1000);
plot(sessionTimes, timelinehWeelPosition(sessionTimeIdx), '*');


photoDiodeTrace = timeline.rawDAQData(:,strcmp(inputNames, 'photoDiode'));
timelineClickTrace = [0;diff(detrend(timeline.rawDAQData(:,strcmp(inputNames, 'audioOut'))))];
for i = 1:length(sessionTimeIdx)
    plotIdx = (sessionTimeIdx(i)-500):(sessionTimeIdx(i)+500);
    plot(timelineTime(plotIdx), photoDiodeTrace(plotIdx)-photoDiodeTrace(plotIdx(1)) + timelinehWeelPosition(plotIdx(1)), 'c');
    plot(timelineTime(plotIdx), timelineClickTrace(plotIdx)-timelineClickTrace(plotIdx(1)) + timelinehWeelPosition(plotIdx(1)), 'm');
end
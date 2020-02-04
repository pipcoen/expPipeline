classdef spatialAnalysis < matlab.mixin.Copyable
    %% spatialAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    % that can be used to plot various aspects of animal behavior. NOTE: This function is designed to operate on the output of convertExpFiles and
    % most methods are specific to multisensoty spatial integration plots.
    %
    % Inputs(default values)
    % subjects(Must be specified)------A cell array of subjects collect parameter and block files for
    % expDate('last')------------------A cell array of dates, one for all subjects or one for each subject
    % combineMice(0)-------------------A tag to indicate whether combine across mice or retain individuals
    %                      0 individual mice, but data is combined across sessions, modal parameter set only
    %                      1 means all mice combined into a single uber-mouse, modal parameter set only
    %                      2 as with 1, but all sessions retained, regardless of conditionParameters used
    % extraTag('eph')-----------------A tag to indicate whether you want to load extra data
    %                      'eph' will load the ephys data, if available, for each session date
    %                      'raw' will load the raw data, (rarely needed wheel times from signals etc.)
    % expDef('multiSpaceWorld')-------Specify the particular expDef to load
    
    properties (Access=public)
        blks;                  %Block files loaded for each subject (one cell per subject)
        glmFit;                %GLM class for post-fitting of data in blks
        hand;                  %Handles to current axis/figure being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate, combineMice, extraTag, expDef)
            % Initialize fields with default values if no vaules are provided. Then called changeMouse function to get requested data.
            if ~exist('subjects', 'var'); error('Must specify which subject to load data from'); end
            if ~exist('expDate', 'var'); expDate = {'last'}; end
            if ~exist('combineMice', 'var'); combineMice = 0; end
            if ~exist('extraTag', 'var'); extraTag = 'eph'; end
            if ~exist('expDef', 'var'); expDef = 'multiSpaceWorld'; end
            
            if ~iscell(expDate); expDate = {expDate}; end
            if ~iscell(subjects); subjects = {subjects}; end
            if any(strcmp(subjects, 'all'))
                subjects = [arrayfun(@(x)['PC0' num2str(x)], 11:99,'uni', 0),'DJ006', 'DJ007','DJ008','DJ010']; 
            end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            expDate = expDate(:); subjects = subjects(:);
            expDate = arrayfun(@(x,y) prc.keyDates(x,y), subjects(:), expDate(:), 'uni', 0);
            subjects = subjects(~cellfun(@isempty, expDate));
            expDate = expDate(~cellfun(@isempty, expDate));
            obj = changeMouse(obj, subjects, expDate, combineMice, expDef, extraTag);
        end
        
        function obj = changeMouse(obj, subjects, expDate, combineMice, extraTag, expDef)
            if ~exist('combineMice', 'var'); combineMice = 0; end
            if ~exist('extraTag', 'var'); extraTag = 'none'; end
            if ~exist('expDef', 'var'); expDef = 'multiSpaceWorld'; end
            obj.blks  = cellfun(@(x,y) prc.getDataFromDates(x, y, extraTag, expDef), subjects(:), expDate(:), 'uni', 0);
            obj.blks = vertcat(obj.blks{:});
            
            mouseList = unique({obj.blks.subject})';
            if combineMice>0; mouseList = {cell2mat(mouseList')}; end
            if combineMice==-2 && length(mouseList)~=1; error('Request one mouse if you want it split into separate days'); end
            if combineMice==-2
                mouseList = arrayfun(@(x) [mouseList{1} ':' num2str(x{1})], {obj.blks.expDate}', 'uni',0);
                [obj.blks.subject] = deal(mouseList{:});
            end
            
            
            retainIdx = ones(length(obj.blks),1)>0;
            subjectRef = zeros(length(obj.blks),1);
            %If there are multiple parameter sets in the requested date range for a mouse, use only the most common parameter set.
            for i = 1:length(mouseList)
                mouseIdx = cellfun(@(x) contains(mouseList{i}, x), {obj.blks.subject}');
                subjectRef(mouseIdx) = i;
                mouseConditions = {obj.blks(mouseIdx).conditionParametersAV}';
                [conditionSets, ~, setIdx] = unique(cellfun(@(x,y) [num2str(x(:)'), y], mouseConditions, {obj.blks(mouseIdx).expType}','uni',0));
                if length(conditionSets)>1 && any(combineMice==[1 0])
                    fprintf('WARNING: Several parameter sets in date range for %s. Using mode\n', mouseList{i});
                    retainIdx(mouseIdx) = retainIdx(mouseIdx).*(setIdx == mode(setIdx));
                end
            end
            obj.blks = obj.blks(retainIdx);
            subjectRef = subjectRef(retainIdx);
            
            obj.blks = cell2mat(arrayfun(@(x) prc.combineBlocks(obj.blks(subjectRef==x)), 1:size(mouseList,1), 'uni', 0)');
            obj.hand.axes = [];
            obj.hand.figure = [];
        end
        
        function viewGLMFit(obj, modelString, cvFolds)
            if ~exist('modelString', 'var'); modelString = 'SimpLog'; end
            if ~exist('cvFolds', 'var'); cvFolds = 0; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.blks);
            axesOpt.btlrMargins = [80 100 80 40];
            axesOpt.gapBetweenAxes = [100 60];
            axesOpt.numOfRows = 2;
            axesOpt.figureHWRatio = 1.1;
            obj.glmFit = cell(length(obj.blks),1);
            for i  = 1:length(obj.blks)
                if contains(lower(modelString), 'nest')
                    normBlock = spatialAnalysis.getBlockType(obj.blks(i),'norm',0);
                    normBlock = prc.filtBlock(normBlock, normBlock.timeOutsBeforeResponse == 0);
                    normBlock.tri.outcome.responseMade(normBlock.responseTime>1.5) = 0;
                    if ~contains(lower(modelString), 'plot'); obj.glmFit{i} = fit.GLMmultiNest(normBlock, modelString); end
                else
                    normBlock = spatialAnalysis.getBlockType(obj.blks(i),'norm');
                    normBlock = prc.filtBlock(normBlock,~isinf(normBlock.tri.stim.audInitialAzimuth));
                    % galvoIdx = ismember(normBlock.galvoPosition, [1.8, -4;3,-4;3,-3], 'rows');
                    % galvoIdx = ismember(normBlock.galvoPosition, [4.2, -2;4.2,-3], 'rows');
                    % galvoIdx = ismember(normBlock.galvoPosition, [0.6, 2; 0.6, 3; 1.8,2], 'rows');
                    % normBlock = prc.filtBlock(normBlock, normBlock.laserType==1 & galvoIdx);
                    if ~contains(lower(modelString), 'plot'); obj.glmFit{i} = fit.GLMmulti(normBlock, modelString); end
                end
                obj.hand.axes = plt.getAxes(axesOpt, i);
                if ~contains(lower(modelString), 'plot'); obj.glmFit{i}.fit; end
                obj.glmFit{i}.plotFit;
                hold on; box off;
                plt.dataWithErrorBars(normBlock,0);
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
                if length(obj.blks) == 1; set(gcf, 'position', get(gcf, 'position').*[1 0.9 1 1.15]); end
                if cvFolds; obj.glmFit{i}.fitCV(cvFolds); end
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
            plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
            obj.hand.figure = [];
        end
        
        
        function getGLMPerturbations(obj, modelString)
            if ~exist('modelString', 'var'); modelString = 'ReducedLogCNSplitDelta'; end
            
            figure;
            initBlock = prc.filtBlock(obj.blks{1}, obj.blks{1}.galvoPosition(:,2)~=4.5);
            initBlock = prc.filtBlock(initBlock, ~ismember(abs(initBlock.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | initBlock.laserType==0);
            initBlock = prc.filtBlock(initBlock, initBlock.timeOutsBeforeResponse==0);
            
            normBlock = spatialAnalysis.getBlockType(initBlock,'norm',1);
            laserBlock = spatialAnalysis.getBlockType(initBlock,'las',1);
            
            obj.glmFit = fit.GLMmulti(normBlock);
            obj.glmFit.GLMMultiModels(modelString);
            obj.glmFit.fit;
            [galvoblks, scanPlot.gridXY] = prc.makeGrid(laserBlock, laserBlock, [], 'galvouni', 3);
            scanPlot.data = nan(size(galvoblks));
            numTrials = zeros(size(galvoblks));
            for j = find(~cellfun(@isempty, galvoblks))'
                if length(galvoblks{j}.responseMade)<500; continue; end
                subGLM = fit.GLMmulti(galvoblks{j});
                subGLM.GLMMultiModels(modelString);
                subGLM.prmInit = obj.glmFit.prmFits;
                subGLM.blockData.freeP = [];
                subGLM.fit;
                minLL = subGLM.calculateLogLik(subGLM.prmFits);
                subGLM.blockData.freeP = 1:length(subGLM.prmFits);%(obj.glmFit.prmFits);
                subGLM.fit;
                maxLL = subGLM.calculateLogLik(subGLM.prmFits);
                
                [xIdx, yIdx] = ind2sub(size(galvoblks), j);
                scanPlot.data(xIdx, yIdx) = maxLL-minLL;
                prmChange(xIdx, yIdx,:) = subGLM.prmFits;
                numTrials(xIdx, yIdx, 1) = length(galvoblks{j}.responseMade);
            end
            %%
            axesOpt = struct('totalNumOfAxes', length(subGLM.prmFits), 'btlrMargins', [50 100 10 10], 'gapBetweenAxes', [40 0], ...
                'numOfRows', 2, 'figureHWRatio', 0.8, 'figureSize', 400);
            for j  = 1:length(subGLM.prmFits)
                hold on; box off;
                scanPlot.data = prmChange(:,:,j);
                scanPlot.title = subGLM.prmLabels{j};
                obj.hand.axes = plt.getAxes(axesOpt, j);
                scanPlot.colorBarLimits = max(abs(scanPlot.data(:)))*[-1 1];
                plt.scanningBrainEffects(scanPlot);
            end
            %             for i  = 1:length(obj.blks)
            %                 hold on; box off;
            %                 obj.hand.axes = plt.getAxes(axesOpt, i);
            %                 scanPlot.colorBarLimits = [-0.1 0.1];
            %                 plt.scanningBrainEffects(scanPlot);
            %             end
        end
        
        function viewDataWithoutFits(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.blks);
            axesOpt.btlrMargins = [80 100 80 40];
            axesOpt.gapBetweenAxes = [100 60];
            axesOpt.numOfRows = ceil(length(obj.blks)/5);
            axesOpt.figureHWRatio = 1.1;
            for i  = 1:length(obj.blks)
                blk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
                plotOpt.Marker = '.'; plotOpt.MarkerSize = 20; plotOpt.lineStyle = '-';
                obj.hand.axes = plt.getAxes(axesOpt, i);
                audValues = unique(blk.exp.conditionParametersAV{1}(:,1));
                visValues = unique(blk.exp.conditionParametersAV{1}(:,2));                
                switch plotType(1:3)
                    case 'rea'
                        gridData = prc.makeGrid(blk, round(blk.tri.outcome.timeToWheelMove*1e3), @median, 1);
                        plt.gridSplitByRows(gridData, visValues*100, audValues, plotOpt);
                    case 'res'
                        gridData = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @mean, 1);
                        plt.gridSplitByRows(gridData, visValues*100, audValues, plotOpt);
                        ylim([0 1]);
                        xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                        yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
                end
                box off;
                title(cell2mat(unique(blk.exp.subject)'));
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
            plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
        end
        
        function scatterData(obj,plotType)
            if ~exist('plotType', 'var'); plotType = 'rea'; end
            figure;
            switch plotType(1:3)
                case 'rea'
                    scatterData = zeros(length(obj.blks),2);
                    for i  = 1:length(obj.blks)
                        nBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm');
                        timeGrid = prc.makeGrid(nBlk, round(nBlk.tri.outcome.timeToWheelMove*1e3), @median, 'abscondition');
                        trialType = nBlk.tri.trialClass;
                        trialType = trialType.blank*0+trialType.auditory + 2*trialType.visual+3*trialType.coherent+4*trialType.conflict;
                        trialGrid = prc.makeGrid(nBlk, trialType, @mean, 'abscondition');
                        tkIdx = trialGrid.*fliplr(trialGrid)==12;
                        scatterData(i,:) = [mean(timeGrid(tkIdx & trialGrid==3)), mean(timeGrid(tkIdx & fliplr(trialGrid)==3))];
                        scatterData(i,:) = scatterData(i,:) - mean(timeGrid(trialGrid==2 | trialGrid==2));
                    end
                    limVal = 60;
                    xLimits = ([-limVal limVal/2]);
                    yLimits = ([-limVal limVal/2]);
                    xlabel2Use = '\fontsize{20} Coherent reaction time (ms)';
                    ylabel2Use = '\fontsize{20} Conflict reaction time (ms)';
                case 'coh'
                    scatterData = zeros(length(obj.blks),2);
                    for i  = 1:length(obj.blks)
                        nBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm');
                        nBlk = prc.filtBlock(nBlk, ~nBlk.tri.trialClass.conflict & ~nBlk.tri.trialClass.blank);
                        perfGrid = prc.makeGrid(nBlk, round(nBlk.tri.outcome.feedbackGiven), @mean, 'abscondition');
                        trialType = nBlk.tri.trialClass;
                        trialType = trialType.blank*0+trialType.auditory + 2*trialType.visual+3*trialType.coherent+4*trialType.conflict;
                        trialGrid = prc.makeGrid(nBlk, trialType, @mean, 'abscondition');
                        columns2use = find(sum(trialGrid == 2 | trialGrid==3)==2);
                        maxUnimodalPer = max([columns2use*0+perfGrid(trialGrid==1); perfGrid(2,columns2use)]);
                        multimodalPer = perfGrid(3,columns2use);
                        scatterData(i,:) = round([mean(maxUnimodalPer) mean(multimodalPer)]*100);
                    end
                    xLimits = ([50 100]);
                    yLimits = ([50 100]);
                    ylabel2Use = '\fontsize{20} Coherent perf (%)';
                    xlabel2Use = '\fontsize{20} Best unimodal perf (%)';
            end
            scatter(scatterData(:,1), scatterData(:,2), 'k', 'markerfacecolor', 'k');
            [~, pVal] = ttest(scatterData(:,1), scatterData(:,2));
            disp(pVal);
            axis equal; axis square;
            xlim(xLimits);
            ylim(yLimits);
            hold on
            plot(xLimits, yLimits, '--k')
            xlabel(xlabel2Use);
            ylabel(ylabel2Use);
        end
        
        
        function viewJitterPlot(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'coh'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.blks);
            axesOpt.btlrMargins = [100 100 100 20];
            axesOpt.gapBetweenAxes = [100 80];
            axesOpt.figureSize = 500;
            for i  = 1:length(obj.blks)
                axesOpt.idx = i;
                switch lower(plotType)
                    case {'coh'; 'con'}
                        plotOpt = prc.getMultiTriplets(normBlock, strcmpi(plotType, 'coh'));
                        plotOpt.figureSize = 400;
                        if strcmpi(plotType, 'coh')
                            plotOpt.mainTitle = '\fontsize{20} Fraction correct by condition';
                            plotOpt.mainYLabel = '\fontsize{20} Fraction correct';
                            plotOpt.yLimits = [0 1.2];
                        else
                            plotOpt.mainTitle = '\fontsize{20} Fraction of unisensory choices in conflict';
                            plotOpt.mainYLabel = '\fontsize{20} Fraction choices in unisensory direction';
                            plotOpt.yLimits = [0 1.2];
                        end
                    case 'whe'
                        normBlock = spatialAnalysis.getBlockType(obj.blks(i),'norm',0);
                        plotOpt = prc.getCoherentConflictPairs(normBlock);
                        plotOpt.mainXLabel = '\fontsize{20} Condition';
                        plotOpt.mainTitle = '\fontsize{20} Difference in reation times for coherent/conflict conditions (ms)';
                        plotOpt.mainYLabel = '\fontsize{20} Time to wheel movement (ms)';
                end
                obj.hand.axes = plt.getAxes(axesOpt);
                plt.jitter(plotOpt.yData, plotOpt); grid('on');
                title(sprintf('%s: n = %d', obj.subjects{i}, obj.blks(i).nSessions));
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel(plotOpt.mainTitle, 't', mainAxes);
            plt.suplabel(plotOpt.mainYLabel, 'y', mainAxes);
            plt.suplabel(plotOpt.mainXLabel, 'x', mainAxes);
        end
        
        
        function runMouseReplicate(obj, subjectNames, funString)
            for i = 1:length(obj.blks)
                replicatedObj = copy(obj);
                replicatedObj.subjects = subjectNames;
                replicatedObj.blks = repmat(replicatedObj.blks(i), length(replicatedObj.subjects),1);
                replicatedObj.expDate = [];
                eval(['replicatedObj.' funString ';']);
                plt.suplabel(obj.subjects{i}, 't');
            end
        end
    end
    
    methods (Static)
        function filteredBlock = getBlockType(block, tag, removeTimeouts)
            if ~exist('tag', 'var'); error('Must specificy tag'); end
            if ~exist('removeTimeouts', 'var'); removeTimeouts = 1; end
            outcome = block.tri.outcome;
            inact = block.tri.inactivation;
            if removeTimeouts; timeOutFilter =  outcome.responseMade~=0; else, timeOutFilter =  outcome.responseMade*0+1; end
            laserOff = prc.filtBlock(block, (inact.laserType==0 | isnan(inact.laserType)) & timeOutFilter & outcome.validTrial, 'tri');
            laserOn = prc.filtBlock(block, (inact.laserType~=0 & ~isnan(inact.laserType)) & timeOutFilter & outcome.validTrial, 'tri');
            switch lower(tag(1:3))
                case 'nor'; filteredBlock = laserOff;
                case 'las'; filteredBlock = laserOn;
            end
        end
        
        function block = removePoorAuditoryDays(block)
            audblks = prc.filtBlock(block, block.trialType==1 & block.laserType == 0 & block.tri.outcome.responseMade~=0);
            sessList = unique(audblks.sessionIdx);
            audPerfomance = arrayfun(@(x) mean(audblks.feedback(audblks.sessionIdx == x)==1), sessList);
            audTrials = arrayfun(@(x) length(audblks.feedback(audblks.sessionIdx == x)==1), sessList);
            block = prc.filtBlock(block, ismember(block.sessionIdx, sessList(audPerfomance>0.7 & audTrials > 50)));
        end
        
        function alterFigure(currentFigureHandle, ~)
            obj = get(currentFigureHandle.Number, 'UserData');
            switch lower(get(currentFigureHandle, 'Tag'))
                case 'boxres'
                    obj.viewBoxPlots('num', currentFigureHandle);
                case 'boxnum'
                    obj.viewBoxPlots('res', currentFigureHandle);
            end
        end
        
        function alterAxes(currentAxesHandle, ~)
            obj = get(currentAxesHandle.Number, 'UserData');
            switch lower(get(currentAxesHandle.Number, 'Tag'))
                case 'boxres'
                    obj.viewBoxPlots('num', currentAxesHandle);
            end
        end
    end
end
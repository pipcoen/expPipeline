classdef spatialAnalysis < matlab.mixin.Copyable
    %% SpatialAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    %that can be used to plot various aspects of animal behavior. NOTE: This function is designed to operate on the output of convertExpFiles and
    %most methods are specific to multisensoty spatial integration plots.
    
    %INPUTS(default values)
    %subjects(required)------------A cell array of subjects collect parameter and block files for
    %expDate('last')---------------A cell array of dates (case insensitive), one for all subjects or one for each subject
    %        'yyyy-mm-dd'--------------A specific date
    %        'all'---------------------All data
    %        'lastx'-------------------The last x days of data (especially useful during training)
    %        'firstx'------------------The first x days of date
    %        'yest'--------------------The x-1 day, where x is the most recent day
    %        'yyyy-mm-dd:yyyy-mm-dd'---Load dates in this range (including the boundaries)
    %        A dataTag (see prc.keyDates)
    %combineMice(0)----------------A tag to indicate whether to combine multiple mice, or split a single mouse into separate days
    %        0-------------------------Individual mice are returned as individual structures
    %        1-------------------------Mice combined into a single uber-mouse with a concatenated name
    %        -1------------------------Only used with one mouse, and means that the mouse is split into individual sessions
    %modalParams(0)----------------Indicates whether to eliminate days where the pamater set (aud/vis combos and expType) didn't match the mode
    %extraTag('eph')---------------Indicates whether you want to load extra data
    %        'eph'---------------------load the ephys data, if available, for each session date
    %        'raw'---------------------load the raw data, (rarely needed wheel times from signals etc.)
    %expDef('multiSpaceWorld')-----Specify the particular expDef to load
    
    properties (Access=public)
        blks;                  %Block files loaded for each subject (one cell per subject)
        glmFit;                %GLM class for post-fitting of data in blks
        hand;                  %Handles to current axis/figure being used for plotting
    end
    
    %%METHODS (there are other methods contained in separate .m files)
    methods
        function obj = spatialAnalysis(subjects, expDate, combineMode, modalParams, extraTag, expDef)
            %% Central function that loads the requested mouse data into a spatial analysis object
            
            %Initialize fields with default values if no vaules are provided.
            if ~exist('subjects', 'var'); error('Must specify which subject to load data from'); end
            if ~exist('expDate', 'var'); expDate = {'last'}; end
            if ~exist('combineMode', 'var'); combineMode = 0; end
            if ~exist('modalParams', 'var'); modalParams = 0; end
            if ~exist('extraTag', 'var'); extraTag = 'eph'; end
            if ~exist('expDef', 'var'); expDef = 'multiSpaceWorld'; end
            
            %Make sure that all fields are cells. If not, convert to cells. If "all" mice are requested, create cell array of all mice.
            if ~iscell(expDate); expDate = {expDate}; end
            if ~iscell(subjects); subjects = {subjects}; end
            if any(strcmp(subjects, 'all'))
                subjects = [arrayfun(@(x)['PC0' num2str(x)], 11:99,'uni', 0),'DJ006', 'DJ007','DJ008','DJ010'];
            end
            
            %If there is only one date provided, repeat for all subjects.
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            subjects = subjects(:); expDate = expDate(:);  %Make sure subjects and rows are columns.
            
            %If a keyDates tag was used (e.g. "behaviour") instead of a "real" date, "prc.keyDates" will get the corresponding date range. If a tag
            %was not used, then the "expDate" input will not match a data tag in prc.keyDates, and so will be unchanged.
            expDate = arrayfun(@(x,y) prc.keyDates(x,y), subjects(:), expDate(:), 'uni', 0);
            subjects = subjects(~cellfun(@isempty, expDate)); %Removes excess subjects in the 'all' case
            expDate = expDate(~cellfun(@isempty, expDate));   %Removes corresponding dates
            
            %Run 'changeMouse' function, which does all the loading and combining of the data
            obj = changeMouse(obj, subjects, expDate, combineMode, modalParams, expDef, extraTag);
        end
        
        function obj = changeMouse(obj, subjects, expDate, combineMode, modalParams, extraTag, expDef)
            %% This function uses the curated inputs to actually load and combine the data as requested
            
            %INPUTS are defined above, but some defaults be redefined here in case this method is called for an existing object
            if ~exist('combineMode', 'var'); combineMode = 0; end
            if ~exist('extraTag', 'var'); extraTag = 'none'; end
            if ~exist('expDef', 'var'); expDef = 'multiSpaceWorld'; end
            
            %Load the data for the requested subjects/dates using prc.getDataFromDates. Concatenate into one structure array, "blks"
            obj.blks  = cellfun(@(x,y) prc.getDataFromDates(x, y, extraTag, expDef), subjects(:), expDate(:), 'uni', 0);
            obj.blks = vertcat(obj.blks{:}); %NOTE: currently errors if some mice have ephys data and others don't
            
            %Get list of subjects and a reference index. Modify depending on the "combineMode"
            mouseList = unique({obj.blks.subject})';
            [~, subjectRef] = ismember({obj.blks.subject}', mouseList);
            if combineMode==1
                mouseList = {cell2mat(mouseList')};
                subjectRef = subjectRef*0+1;
            elseif combineMode==-1 
                if length(mouseList)~=1; error('Request one mouse if you want it split into separate days'); end
                mouseList = arrayfun(@(x) [mouseList{1} ':' num2str(x{1})], {obj.blks.expDate}', 'uni',0);
                [obj.blks.subject] = deal(mouseList{:});
                subjectRef = (1:length(subjectRef))';
            end
            
            %This loop removes the less common sessions for each mouse if "modalParams==1" and the mouse has a mix of paramters/expTypes
            retainIdx = ones(length(obj.blks),1)>0;
            for i = 1:max(subjectRef)
                mouseConditions = {obj.blks(subjectRef==i).conditionParametersAV}';
                mouseExpTypes = {obj.blks(subjectRef==i).expType}';
                [conditionSets, ~, setIdx] = unique(cellfun(@(x,y) [num2str(x(:)'), y], mouseConditions, mouseExpTypes,'uni',0));
                if length(conditionSets)>1 && modalParams
                    fprintf('WARNING: Several parameter sets in date range for %s. Using mode\n', mouseList{i});
                    retainIdx(subjectRef==i) = retainIdx(subjectRef==i).*(setIdx == mode(setIdx));
                end
            end
            obj.blks = obj.blks(retainIdx);
            subjectRef = subjectRef(retainIdx);
            
            %Run the blk array through "combineBlocks" which converts blks into their final format
            obj.blks = cell2mat(arrayfun(@(x) prc.combineBlocks(obj.blks(subjectRef==x)), 1:size(mouseList,1), 'uni', 0)');
            obj.hand.axes = [];
            obj.hand.figure = [];
        end
               
        function getGLMPerturbations(obj, modelString)
            if ~exist('modelString', 'var'); modelString = 'ReducedLogCN'; end
            
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
        
        function scatterData(obj,plotType)
            if ~exist('plotType', 'var'); plotType = 'rea'; end
            figure;
            switch plotType(1:3)
                case 'rea'
                    scatterData = zeros(length(obj.blks),2);
                    for i  = 1:length(obj.blks)
                        nBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm');
                        timeGrid = prc.makeGrid(nBlk, round((nBlk.tri.outcome.timeToFirstMove-nBlk.tri.timings.stimPeriodStart)*1e3), @nanmean, 'abscondition');
                        trialType = nBlk.tri.trialClass;
                        trialType = trialType.blank*0+trialType.auditory + 2*trialType.visual+3*trialType.coherent+4*trialType.conflict;
                        trialGrid = prc.makeGrid(nBlk, trialType, @mean, 'abscondition');
                        tkIdx = trialGrid.*fliplr(trialGrid)==12;
                        scatterData(i,:) = [mean(timeGrid(tkIdx & trialGrid==3)), mean(timeGrid(tkIdx & fliplr(trialGrid)==3))];
                        scatterData(i,:) = scatterData(i,:) - mean(timeGrid(trialGrid==2 | trialGrid==2));
                    end
                    limVal = 100;
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
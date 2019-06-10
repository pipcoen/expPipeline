classdef spatialAnalysis < matlab.mixin.Copyable
    %% spatialAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    % that can be used to plot various aspects of animal behavior. NOTE: This function is designed to operate on the output of convertExpFiles and
    % most methods are specific to multisensoty spatial integration plots.
    %
    % Inputs(default values) subjects({'PC011';'PC012';'PC013';'PC016'})------A c   ell array of subjects collect parameter and block files for
    % expDate('last')----------------------------------A cell array of dates, one for all subjects or one for each subject
    % combineMice(0)----------------------------A tag to indicate whether you want to combine data across mice or retain individuals.
    
    properties (Access=public)
        subjects;                %Cell array of subject names
        expDate;                 %Recording dates to use for each subject
        blocks;                  %Block files loaded for each subject (one cell per subject)
        glmFit;
        hand;                    %Handles to current axis/figure being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate, combineMice, extraTag, specificExpDef)
            % Initialize fields with default values if no vaules are provided. Then called changeMouse function to get requested data.
            prc.updatePaths;
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~iscell(expDate); expDate = {expDate}; end
            if ~exist('combineMice', 'var'); combineMice = 0; end
            if ~exist('extraTag', 'var'); extraTag = 'none'; end
            if ~exist('specificExpDef', 'var'); specificExpDef = 'multiSpaceWorld'; end
            if ~exist('subjects', 'var') || any(strcmp(subjects, 'all')); subjects = prc.keyDates({'all'}, expDate); end
            if ~iscell(subjects); subjects = {subjects}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            expDate = expDate(:); subjects = subjects(:);
            expDate = arrayfun(@(x,y) prc.keyDates(x,y), subjects(:), expDate(:), 'uni', 0);
            obj = changeMouse(obj, subjects, expDate, combineMice, extraTag, specificExpDef);
        end
        
        
        function viewGLMFit(obj, modelString, cvFolds)
            if ~exist('modelString', 'var'); modelString = 'SimpLog'; end
            if ~exist('cvFolds', 'var'); cvFolds = 0; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [80 100 80 40];
            axesOpt.gapBetweenAxes = [100 60];
            axesOpt.numOfRows = 2;
            axesOpt.figureHWRatio = 1.1;
            obj.glmFit = cell(length(obj.subjects),1);
            for i  = 1:length(obj.subjects)
                if contains(lower(modelString), 'nest')
                    normBlock = spatialAnalysis.getBlockType(obj.blocks{i},'norm',0);
                    normBlock = prc.filtStruct(normBlock, normBlock.timeOutsBeforeResponse == 0);
                    normBlock.responseMade(normBlock.responseTime>1.5) = 0;
                    if ~contains(lower(modelString), 'plot'); obj.glmFit{i} = fit.GLMmultiNest(normBlock, modelString); end
                else
                    normBlock = spatialAnalysis.getBlockType(obj.blocks{i},'norm');
                    % galvoIdx = ismember(normBlock.galvoPosition, [1.8, -4;3,-4;3,-3], 'rows');
                    % galvoIdx = ismember(normBlock.galvoPosition, [4.2, -2;4.2,-3], 'rows');
                    % galvoIdx = ismember(normBlock.galvoPosition, [0.6, 2; 0.6, 3; 1.8,2], 'rows');
                    % normBlock = prc.filtStruct(normBlock, normBlock.laserType==1 & galvoIdx);
                    if ~contains(lower(modelString), 'plot'); obj.glmFit{i} = fit.GLMmulti(normBlock, modelString); end
                end
                obj.hand.axes = plt.getAxes(axesOpt, i);
                if ~contains(lower(modelString), 'plot'); obj.glmFit{i}.fit; end
                obj.glmFit{i}.plotFit;
                hold on; box off;
                plt.dataWithErrorBars(normBlock,0);
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
                if length(obj.subjects) == 1; set(gcf, 'position', get(gcf, 'position').*[1 0.9 1 1.15]); end
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
            initBlock = prc.filtStruct(obj.blocks{1}, obj.blocks{1}.galvoPosition(:,2)~=4.5);
            initBlock = prc.filtStruct(initBlock, ~ismember(abs(initBlock.galvoPosition(:,1)),[0.5; 2; 3.5; 5]) | initBlock.laserType==0);
            initBlock = prc.filtStruct(initBlock, initBlock.timeOutsBeforeResponse==0);
            
            normBlock = spatialAnalysis.getBlockType(initBlock,'norm',1);
            laserBlock = spatialAnalysis.getBlockType(initBlock,'las',1);
            
            obj.glmFit = fit.GLMmulti(normBlock);
            obj.glmFit.GLMMultiModels(modelString);
            obj.glmFit.fit;
            [galvoBlocks, scanPlot.gridXY] = prc.makeGrid(laserBlock, laserBlock, [], 'galvouni', 3);
            scanPlot.data = nan(size(galvoBlocks));
            numTrials = zeros(size(galvoBlocks));
            for j = find(~cellfun(@isempty, galvoBlocks))'
                if length(galvoBlocks{j}.responseMade)<500; continue; end
                subGLM = fit.GLMmulti(galvoBlocks{j});
                subGLM.GLMMultiModels(modelString);
                subGLM.prmInit = obj.glmFit.prmFits;
                subGLM.blockData.freeP = [];
                subGLM.fit;
                minLL = subGLM.calculateLogLik(subGLM.prmFits);
                subGLM.blockData.freeP = 1:length(subGLM.prmFits);%(obj.glmFit.prmFits);
                subGLM.fit;
                maxLL = subGLM.calculateLogLik(subGLM.prmFits);
                
                [xIdx, yIdx] = ind2sub(size(galvoBlocks), j);
                scanPlot.data(xIdx, yIdx) = maxLL-minLL;
                prmChange(xIdx, yIdx,:) = subGLM.prmFits;
                numTrials(xIdx, yIdx, 1) = length(galvoBlocks{j}.responseMade);
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
%             for i  = 1:length(obj.subjects)
%                 hold on; box off;
%                 obj.hand.axes = plt.getAxes(axesOpt, i);
%                 scanPlot.colorBarLimits = [-0.1 0.1];
%                 plt.scanningBrainEffects(scanPlot);
%             end
        end
        
        function viewDataWithoutFits(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [80 100 80 40];
            axesOpt.gapBetweenAxes = [100 60];
            axesOpt.numOfRows = ceil(length(obj.subjects)/5);
            axesOpt.figureHWRatio = 1.1;
            for i  = 1:length(obj.subjects)
                normBlock = spatialAnalysis.getBlockType(obj.blocks{i},'norm',1);
                plotOpt.Marker = '.'; plotOpt.MarkerSize = 20; plotOpt.lineStyle = '-';
                obj.hand.axes = plt.getAxes(axesOpt, i);
                switch plotType(1:3)
                    case 'rea'
                        gridData = prc.makeGrid(normBlock, round(normBlock.timeToWheelMove*1e3), @median, 1);
                        plt.gridSplitByRows(gridData, normBlock.visValues*100, normBlock.audValues, plotOpt);
                    case 'res'
                        gridData = prc.makeGrid(normBlock, normBlock.responseMade==2, @mean, 1);
                        plt.gridSplitByRows(gridData, normBlock.visValues*100, normBlock.audValues, plotOpt);
                        ylim([0 1]);
                        xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                        yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
                end
                box off;
                title(obj.subjects{i});
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
            plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
        end
        
        function scatterReactionTimes(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            allTimes = [];
            for i  = 1:length(obj.subjects)
                normBlock = spatialAnalysis.getBlockType(obj.blocks{i},'norm');
                timeGrid = prc.makeGrid(normBlock, round(normBlock.timeToWheelMove*1e3), @median, 'abscondition');
                trialGrid = prc.makeGrid(normBlock, round(normBlock.trialType), @mean, 'abscondition');
                tkIdx = trialGrid.*fliplr(trialGrid)==12;
                allTimes = [allTimes; [mean(timeGrid(tkIdx & trialGrid==3)), mean(timeGrid(tkIdx & fliplr(trialGrid)==3))]];
%                 allTimes(end,:) = allTimes(end,:)./nanmean(timeGrid(:));
            end
            scatter(allTimes(:,1), allTimes(:,2), 'k', 'markerfacecolor', 'k');
            [~, pVal] = ttest(allTimes(:,1), allTimes(:,2));
            disp(pVal);
            xlim([200 450]);
            ylim([200 450]);
            axis equal; axis square;
            hold on
            plot([200,500], [200,500], '--k')
            xlabel('\fontsize{20} Coherent reaction time (ms)');
            ylabel('\fontsize{20} Conflict reaction time (ms)');
        end
        
        function classifyNeurons(obj)
            clusterSigLevel = kil.findResponsiveCells(obj.blocks{1});
            
        end
        
        function viewJitterPlot(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'coh'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [100 100 100 20];
            axesOpt.gapBetweenAxes = [100 80];
            axesOpt.figureSize = 500;
            for i  = 1:length(obj.subjects)
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
                        normBlock = spatialAnalysis.getBlockType(obj.blocks{i},'norm',0);
                        plotOpt = prc.getCoherentConflictPairs(normBlock);
                        plotOpt.mainXLabel = '\fontsize{20} Condition';
                        plotOpt.mainTitle = '\fontsize{20} Difference in reation times for coherent/conflict conditions (ms)';
                        plotOpt.mainYLabel = '\fontsize{20} Time to wheel movement (ms)';
                end
                obj.hand.axes = plt.getAxes(axesOpt);
                plt.jitter(plotOpt.yData, plotOpt); grid('on');
                title(sprintf('%s: n = %d', obj.subjects{i}, obj.blocks{i}.nSessions));
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel(plotOpt.mainTitle, 't', mainAxes);
            plt.suplabel(plotOpt.mainYLabel, 'y', mainAxes);
            plt.suplabel(plotOpt.mainXLabel, 'x', mainAxes);
        end
        
        
        function obj = changeMouse(obj, subjects, expDate, combineMice, extraTag, specificExpDef)
            if ~exist('combineMice', 'var'); combineMice = 0; end
            if ~exist('extraTag', 'var'); extraTag = []; end
            if ~exist('specificExpDef', 'var'); specificExpDef = 'multiSpaceWorld'; end
            dataType = [extraTag 'bloprm'];
            obj.blocks  = arrayfun(@(x,y) prc.getFilesFromDates(x, y{1}, dataType, specificExpDef), subjects(:), expDate(:), 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.expDate = expDate;
            
            mouseList = unique({obj.blocks.subject})';
            if combineMice==1; mouseList = {cell2mat(mouseList')}; end   
            if combineMice==-2 && length(mouseList)~=1; error('Request one mouse if you want it split into separate days'); end
            if combineMice==-2
                mouseList = arrayfun(@(x) [mouseList{1} ':' num2str(x{1})], {obj.blocks.expDate}', 'uni',0);
                obj.subjects = mouseList;
                obj.expDate = {obj.blocks.expDate}';
                [obj.blocks.subject] = deal(mouseList{:});
            end
                
            
            retainIdx = ones(length(obj.blocks),1)>0;
            %If there are multiple parameter sets in the requested date range for a mouse, use only the most common parameter set.
            for i = 1:length(mouseList)
                mouseIdx = cellfun(@(x) contains(mouseList{i}, x), {obj.blocks.subject}');
                mouseParams = [obj.blocks(mouseIdx).params];
                mouseConditions = {obj.blocks(mouseIdx).uniqueConditions}';
                mouseConditions = cellfun(@(x) [x(:,1)>0 x(:,2:end)], mouseConditions, 'uni', 0); %ignoring different aud amplitudes for now
                [conditionSets, ~, setIdx] = unique(cellfun(@(x,y) num2str([x(:)', y]), mouseConditions, {mouseParams.laserSession}','uni',0));
                if length(conditionSets)>1 && combineMice>=0
                    fprintf('WARNING: Several parameter sets in date range for %s. Using mode\n', mouseList{i});
                    retainIdx(mouseIdx) = retainIdx(mouseIdx).*(setIdx == mode(setIdx));
                end
            end
            obj.blocks = obj.blocks(retainIdx);
            
            obj.subjects = unique({obj.blocks.subject})';
            if combineMice==1; obj.expDate = cell2mat([obj.expDate{ismember(subjects, obj.subjects)}]); obj.subjects = {cell2mat(obj.subjects')}; end          
            [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
            if combineMice==1; subjectIdx = subjectIdx*0+1; end

            idx = 1:length(obj.subjects);
            obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x)), idx, 'uni', 0)';
            obj.hand.axes = [];
            obj.hand.figure = [];
        end
        
        function runMouseReplicate(obj, subjectNames, funString)
            for i = 1:length(obj.subjects)
                replicatedObj = copy(obj);
                replicatedObj.subjects = subjectNames;
                replicatedObj.blocks = repmat(replicatedObj.blocks(i), length(replicatedObj.subjects),1);
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
            if removeTimeouts; timeOutFilter =  block.responseMade~=0; else, timeOutFilter =  block.responseMade*0+1; end
            laserOff = prc.filtStruct(block, (block.laserType==0 | isnan(block.laserType)) & timeOutFilter);
            laserOn = prc.filtStruct(block, (block.laserType~=0 & ~isnan(block.laserType)) & timeOutFilter);            
            switch lower(tag(1:3))
                case 'nor'; filteredBlock = laserOff;
                case 'las'; filteredBlock = laserOn;
            end
        end
        
        function block = removePoorAuditoryDays(block)
            audBlocks = prc.filtStruct(block, block.trialType==1 & block.laserType == 0 & block.responseMade~=0);
            sessList = unique(audBlocks.sessionIdx);
            audPerfomance = arrayfun(@(x) mean(audBlocks.feedback(audBlocks.sessionIdx == x)==1), sessList);
            audTrials = arrayfun(@(x) length(audBlocks.feedback(audBlocks.sessionIdx == x)==1), sessList);
            block = prc.filtStruct(block, ismember(block.sessionIdx, sessList(audPerfomance>0.7 & audTrials > 50)));
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
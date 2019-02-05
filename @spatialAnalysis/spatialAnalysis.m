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
        params;                  %Parameter files loaded for each subject (one cell per subject)
        glmFit;
        axesHandles;             %Handle to current axis being used for plotting
        figureHandles;           %Handle to current figure being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate, combineMice, extraTag)
            % Initialize fields with default values if no vaules are provided. Then called changeMouse function to get requested data.
            prc.updatePaths;
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~iscell(expDate); expDate = {expDate}; end
            if ~exist('combineMice', 'var'); combineMice = 0; end
            if ~exist('extraTag', 'var'); extraTag = 'none'; end
            if ~exist('subjects', 'var') || any(strcmp(subjects, 'all')); subjects = prc.keyDates({'all'}, expDate); end
            if ~iscell(subjects); subjects = {subjects}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            expDate = expDate(:); subjects = subjects(:);
            expDate = arrayfun(@(x,y) prc.keyDates(x,y), subjects(:), expDate(:), 'uni', 0);
            obj = changeMouse(obj, subjects, expDate, combineMice, extraTag);
        end
        
        
        function viewGLMFit(obj, modelString, cvFolds)
            if ~exist('modelString', 'var'); modelString = 'SqrtLogisticSplitNest'; end
            if ~exist('cvFolds', 'var'); cvFolds = 0; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [80 100 80 40];
            axesOpt.gapBetweenAxes = [100 60];
            axesOpt.numOfRows = 2;
            axesOpt.figureHWRatio = 1.1;
            for i  = 1:length(obj.subjects)
                if contains(lower(modelString), 'nest')
                    [normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i}, 0, 5);
                    normBlock = prc.filtStruct(normBlock, normBlock.timeOutsBeforeResponse == 0);
                    normBlock.responseMade(normBlock.responseTime>1.5) = 0;
                    if ~contains(lower(modelString), 'plot'); obj.glmFit{i} = fit.GLMmultiNest(normBlock, modelString); end
                else
                    %                     [~, normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                    [normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                    %                     galvoIdx = ismember(normBlock.galvoPosition, [1.8, -4;3,-4;3,-3], 'rows');
                    %                     galvoIdx = ismember(normBlock.galvoPosition, [4.2, -2;4.2,-3], 'rows');
                    %                     galvoIdx = ismember(normBlock.galvoPosition, [0.6, 2; 0.6, 3; 1.8,2], 'rows');
                    %                     normBlock = prc.filtStruct(normBlock, normBlock.laserType==1 & galvoIdx);
                    if ~contains(lower(modelString), 'plot'); obj.glmFit{i} = fit.GLMmulti(normBlock, modelString); end
                end
                obj.axesHandles = plt.getAxes(axesOpt, i);
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
            obj.figureHandles = [];
        end
        
        function getGLMPerturbations(obj, modelString)
            if ~exist('modelString', 'var'); modelString = 'SimpLogSplitDelta'; end
            axesOpt = struct('totalNumOfAxes', length(obj.subjects), 'btlrMargins', [50 100 10 10], 'gapBetweenAxes', [40 0], ...
                'numOfRows', 2, 'figureHWRatio', 0.8, 'figureSize', 400);
            
            figure;
            [normBlock, laserBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{1}, 1);
            normBlock = prc.filtStruct(normBlock, normBlock.galvoPosition(:,2)~=4.5);
            laserBlock = prc.filtStruct(laserBlock, laserBlock.galvoPosition(:,2)~=4.5 & laserBlock.laserType == 1);
            obj.glmFit = fit.GLMmulti(normBlock);
            obj.glmFit.GLMMultiModels(modelString);
            obj.glmFit.fit;
            [galvoBlocks, scanPlot.gridXY] = prc.makeGrid(laserBlock, laserBlock, [], 'galvouni', 2);
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
                subGLM.blockData.freeP = [1:3 5:8];%(obj.glmFit.prmFits);
                subGLM.fit;
                maxLL = subGLM.calculateLogLik(subGLM.prmFits);
                
                [xIdx, yIdx] = ind2sub(size(galvoBlocks), j);
                scanPlot.data(xIdx, yIdx) = maxLL-minLL;
                numTrials(xIdx, yIdx, 1) = length(galvoBlocks{j}.responseMade);
            end
            for i  = 1:length(obj.subjects)
                hold on; box off;
                obj.axesHandles = plt.getAxes(axesOpt, i);
                scanPlot.colorBarLimits = [-0.1 0.1];
                plt.scanningBrainEffects(scanPlot);
            end
        end
        
        function viewDataWithoutFits(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [80 100 80 40];
            axesOpt.gapBetweenAxes = [100 60];
            axesOpt.numOfRows = 2;
            axesOpt.figureHWRatio = 1.1;
            for i  = 1:length(obj.subjects)
                %                 goodBlock = spatialAnalysis.removePoorAuditoryDays(obj.blocks{i});
                normBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                plotOpt.Marker = '.'; plotOpt.MarkerSize = 20; plotOpt.lineStyle = '-';
                obj.axesHandles = plt.getAxes(axesOpt, i);
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
                normBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                timeGrid = prc.makeGrid(normBlock, round(normBlock.timeToWheelMove*1e3), @median, 'abscondition');
                trialGrid = prc.makeGrid(normBlock, round(normBlock.trialType), @mean, 'abscondition');
                tkIdx = trialGrid.*fliplr(trialGrid)==12;
                allTimes = [allTimes; [mean(timeGrid(tkIdx & trialGrid==3)), mean(timeGrid(tkIdx & fliplr(trialGrid)==3))]];
            end
            scatter(allTimes(:,1), allTimes(:,2), 'k', 'markerfacecolor', 'k');
            xlim([200 450]);
            ylim([200 450]);
            axis equal; axis square;
            hold on
            plot([200,500], [200,500], '--k')
            xlabel('\fontsize{20} Coherent reaction time (ms)');
            ylabel('\fontsize{20} Conflict reaction time (ms)');
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
                        normBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                        plotOpt = prc.getCoherentConflictPairs(normBlock);
                        plotOpt.mainXLabel = '\fontsize{20} Condition';
                        plotOpt.mainTitle = '\fontsize{20} Difference in reation times for coherent/conflict conditions (ms)';
                        plotOpt.mainYLabel = '\fontsize{20} Time to wheel movement (ms)';
                end
                obj.axesHandles = plt.getAxes(axesOpt);
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
        
        
        function obj = changeMouse(obj, subjects, expDate, combineMice, extraTag)
            if ~exist('combineMice', 'var'); combineMice = 0; end
            dataType = [extraTag 'bloprm'];
            [obj.blocks, obj.params]  = arrayfun(@(x,y) prc.getFilesFromDates(x, y{1}, dataType), subjects(:), expDate(:), 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
            obj.expDate = expDate;
            
            mouseList = unique({obj.blocks.subject})';
            if combineMice; mouseList = {cell2mat(mouseList')}; end       
            
            retainIdx = ones(length(obj.params),1)>0;
            %If there are multiple parameter sets in the requested date range for a mouse, use only the most common parameter set.
            for i = 1:length(mouseList)
                mouseIdx = cellfun(@(x) contains(mouseList{i}, x), {obj.blocks.subject}');
                [conditionSets, ~, setIdx] = unique(cellfun(@(x,y) num2str([x(:)', unique(y(:))]),{obj.blocks(mouseIdx).uniqueConditions}', {obj.blocks(mouseIdx).laserSession}','uni',0));
                if length(conditionSets)>1
                    fprintf('WARNING: Several parameter sets in date range for %s. Using mode\n', mouseList{i});
                    retainIdx(mouseIdx) = retainIdx(mouseIdx).*(setIdx == mode(setIdx));
                end
            end
            obj.blocks = obj.blocks(retainIdx);
            obj.params = obj.params(retainIdx);
            
            obj.subjects = unique({obj.blocks.subject})';
            if combineMice; obj.expDate = cell2mat([obj.expDate{ismember(subjects, obj.subjects)}]); obj.subjects = {cell2mat(obj.subjects')}; end          
            [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
            if combineMice; subjectIdx = subjectIdx*0+1; end

            idx = 1:length(obj.subjects);
            obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x)), idx, 'uni', 0)';
            obj.params = arrayfun(@(x) prc.combineParams(obj.params(subjectIdx==x)), idx, 'uni', 0)';
            obj.axesHandles = [];
            obj.figureHandles = [];
        end
        
        function runMouseReplicate(obj, subjectNames, funString)
            for i = 1:length(obj.subjects)
                replicatedObj = copy(obj);
                replicatedObj.subjects = subjectNames;
                replicatedObj.blocks = repmat(replicatedObj.blocks(i), length(replicatedObj.subjects),1);
                replicatedObj.params = repmat(replicatedObj.params(i), length(replicatedObj.subjects),1);
                replicatedObj.expDate = [];
                eval(['replicatedObj.' funString ';']);
                plt.suplabel(obj.subjects{i}, 't');
            end
        end
    end
    
    methods (Static)
        function [normBlock, laserBlock] = getMaxNumberOfTrials(block, inactivation, excludedResponses)
            if ~exist('inactivation', 'var'); inactivation = 0; end
            if ~exist('excludedResponses', 'var'); excludedResponses = 0; end
            if inactivation == 2
                normBlock = prc.filtStruct(block, block.laserSession~=0 & block.responseMade~=excludedResponses);
            elseif mode(block.laserSession>0)==1 || inactivation
                normBlock = prc.filtStruct(block, block.laserSession~=0 & block.laserType==0 & block.responseMade~=excludedResponses);
                laserBlock = prc.filtStruct(block, block.laserSession~=0 & block.laserType~=0 & block.responseMade~=excludedResponses);
            else, normBlock = prc.filtStruct(block, block.laserSession==0 & block.responseMade~=excludedResponses);
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
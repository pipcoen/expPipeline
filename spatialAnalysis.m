classdef spatialAnalysis
    %% spatialAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    % that can be used to plot various aspects of animal behavior. NOTE: This function is designed to operate on the output of convertExpFiles and
    % most methods are specific to multisensoty spatial integration plots.
    %
    % Inputs(default values) subjects({'PC011';'PC012';'PC013';'PC016'})------A cell array of subjects collect parameter and block files for
    % expDate('last')----------------------------------A cell array of dates, one for all subjects or one for each subject
    % processingTag('none')----------------------------A tag to indicate whether you want to combine data across mice or retain individuals.
    
    properties (Access=public)
        subjects;                %Cell array of subject names
        expDate;                 %Recording dates to use for each subject
        blocks;                  %Block files loaded for each subject (one cell per subject)
        params;                  %Parameter files loaded for each subject (one cell per subject)
        splits;                  %Block splits blocks that are commonly used
        currentAxes;             %Handle to current axis being used for plotting
        currentFigure;           %Handle to current figure being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate)
            % Initialize fields with default values if no vaules are provided.
            if ~exist('subjects', 'var') || isempty(subjects); subjects = {'PC011';'PC012';'PC013';'PC015';'PC010';'PC017'}; end
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~iscell(subjects); subjects = {subjects}; end
            if ~iscell(expDate); expDate = {expDate}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            obj = changeMouse(obj, subjects, expDate);
        end
        
        function viewBoxPlots(obj, plotType)
            figure;
            if ~exist('plotType', 'var'); plotType = 'res'; end
            maxGrid = max(cell2mat(cellfun(@(x) [length(x.audValues) length(x.visValues)], obj.blocks, 'uni', 0)), [], 1);
            colorYTick = [0 1]; colorMap = plt.redblue(64); colorDirection = 'normal'; axisLimits = [0 1]; colorLabel = 'Fraction of right turns';
                        
            for i  = 1:length(obj.subjects)
                if contains(plotType(1:3), {'las'; 'ina'}) || mode(obj.blocks{i}.laserSession>0)==1
                    laserOff = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserSession~=0 & obj.blocks{i}.laserType==0);
                else, laserOff = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserSession==0);
                end
                
                switch lower(plotType(1:3))
                    case 'res'
                        plotData = prc.makeGrid(laserOff, obj.splits.maxNoLaser{i}.response==2, @mean, 1);
                        trialNumber = length(laserOff.response);
                    case {'svd'; 'mod'}
                        if ~isempty(obj.expDate)
                            obj.runMouseReplicate({'Original'; 'Model'; 'Difference'}, ['viewBoxPlots(''' plotType ''')']); 
                            return;
                        end
                        numRightTurns = prc.makeGrid(laserOff, laserOff.response==2, @sum, 1);
                        numTrials = prc.makeGrid(laserOff, laserOff.response~=-1, @sum, 1);
                        responseData = numRightTurns./numTrials;
                        if strcmp(plotType, 'svd')
                            obj.blocks{i}.visValues = obj.blocks{i}.visValues(~any(isnan(responseData)));
                            responseData = responseData(:,~any(isnan(responseData),1));
                            [U,S,V] = svd(responseData, 'econ');
                            modelData = U(:,1)*V(:,1)'*S(1,1);
                        elseif strcmp(plotType, 'mod')
                            audIdx = obj.blocks{1}.grids.visValues==0;
                            visIdx = obj.blocks{1}.grids.audValues==0;
                            modelData = sqrt((numRightTurns(audIdx)*numRightTurns(visIdx)')./(numTrials(audIdx)*numTrials(visIdx)'));
                            modelData(isnan(responseData)) = nan;
                        end
                        
                        axisLimits = [-0.25 1.25];
                        trialNumber = length(laserOff.response);
                        if i == 1; plotData = responseData; end
                        if i == 2; plotData = modelData; end
                        if i == 3; plotData = responseData - modelData; axisLimits = [-0.25 0.25]; end
                end
                plotIdx = ~all(isnan(plotData));
                plt.getAxes(i, length(obj.subjects), [], maxGrid(2)/(1.3*maxGrid(1)), [100 80 60 100], [100 40]);
                plotData = plotData(:,plotIdx);
                imsc(plotData, axisLimits, colorMap, 'k');
                [xPnts, yPnts] = meshgrid(1:size(plotData,2), 1:size(plotData,1));
                arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center'), xPnts, yPnts, plotData) 
                daspect([1 1 1]); axis xy;
                title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, trialNumber, obj.blocks{i}.nSessions));
                set(gca, 'xTick', 1:size(plotData,2), 'xTickLabel', obj.blocks{i}.visValues(plotIdx)*100, 'fontsize', 14)
                set(gca, 'yTick', 1:size(plotData,1), 'yTickLabel', obj.blocks{i}.audValues, 'fontsize', 14, 'TickLength', [0, 0])
                ylabel(obj.blocks{i}.audType);
                xlabel('Visual Contast');
                box off;
            end
            currentAxisPotision = get(gca, 'position');
            colorBar = colorbar;
            set(colorBar,'Ticks', get(colorBar, 'Limits'), 'TickLabels', colorYTick, 'YDir', colorDirection);
            set(gca, 'position', currentAxisPotision);
            figureSize = get(gcf, 'position');
            colorLabel = ylabel(colorBar, colorLabel);
            set(colorLabel, 'position', [1 mean(get(colorBar, 'Limits')) 0], 'FontSize', 14)
            set(colorBar, 'position', [1-75/figureSize(3), 0.2, 30/figureSize(3), 0.6])
        end
        
        function viewGLMFit(obj, modelString)
            if ~exist('modelString', 'var'); modelString = 'C-subset-multiAud2'; end
            figure;
            for i  = 1:length(obj.subjects)
                [obj.currentAxes, obj.currentFigure] = plt.getAxes(i, length(obj.subjects), [], [], [60 80 60 20], [200 40]);
                glm = fit.GLMmulti(obj.blocks{i});
                glm = glm.setModel(modelString);
                glm = glm.fit;
                glm.plotFit;              
            end
        end
        
        function viewPsychoCurves(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            for i  = 1:length(obj.subjects)
                audValues = obj.blocks{i}.audValues;
                colorChoices = plt.selectRedBlueColors(audValues);
                visGrid = obj.blocks{i}.grids.visValues;
                audGrid = obj.blocks{i}.grids.audValues;
                switch plotType(1:3)
                    case 'res'
                        obj.currentAxes = plt.getAxes(i, length(obj.subjects), [], [], [60 80 60 20], [200 40]); 
                        hold on; box off;
                        numTrials = prc.makeGrid(obj.blocks{i}, obj.blocks{i}.laserPower==0, @length, 1);
                        numRightTurns = prc.makeGrid(obj.blocks{i}, obj.blocks{i}.laserPower==0 & obj.blocks{i}.response==2, @sum, 1);
                        for j = 1:length(audValues)
                            idx = audGrid==audValues(j) & numTrials>0;
                            StimLevelsFineGrain = min(visGrid(:)):(range(visGrid(:))/1000):max(visGrid(:));
                            [paramsValues, fittingFunction] = fit.psychoCurve(visGrid(idx)', numRightTurns(idx), numTrials(idx));
                            plot(StimLevelsFineGrain, fittingFunction(paramsValues, StimLevelsFineGrain),'LineWidth',1.1, 'color', colorChoices(j,:));
                            plot(visGrid(idx), numRightTurns(idx)./numTrials(idx), '.', 'color', colorChoices(j,:));
                        end
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, trialEst, laserOff.nSessions));
                        set(gca, 'XTick', -45:15:45);
                end
            end
            figureSize = get(gcf, 'position');
            mainAxes = [60./figureSize(3:4) 1-2*(50./figureSize(3:4))];
            plt.suplabel('\fontsize{16} Auditory on {\color{red}right, \color{black}centre, or \color{blue}left}', 't', mainAxes);
            plt.suplabel(['\fontsize{16} ' yAxLabel], 'y', mainAxes);
            plt.suplabel('\fontsize{16} Visual Contrast', 'x', mainAxes);
        end    
        
        function obj = changeMouse(obj, subjects, expDate)
            %Get block and parameter files for the requested dates.
            [obj.blocks, obj.params]  = cellfun(@(x) prc.getFilesFromDates(x, expDate, 'bloprm'), subjects, 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
            obj.expDate = expDate;
            
            retainIdx = ones(length(obj.params),1)>0;
            %If there are multiple parameter sets in the requested date range for a mouse, use only the most common parameter set.
            for i = 1:length(subjects)
                mouseIdx = strcmp({obj.blocks.subject}', subjects{i});
                [conditionSets, ~, setIdx] = unique(cellfun(@(x) num2str(x(:)'),{obj.blocks(mouseIdx).uniqueConditions}','uni',0));
                if length(conditionSets)>1
                    fprintf('WARNING: Several parameter sets in date range for %s. Using mode\n', subjects{i});
                    retainIdx(mouseIdx) = retainIdx(mouseIdx).*(setIdx == mode(setIdx));
                end
            end
            obj.blocks = obj.blocks(retainIdx);
            obj.params = obj.params(retainIdx);
            obj.subjects = unique({obj.blocks.subject}');
            
            [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
            idx = 1:length(subjects);
            obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x)), idx, 'uni', 0)';
            obj.params = arrayfun(@(x) prc.combineParams(obj.params(subjectIdx==x)), idx, 'uni', 0)';
            
            obj.currentAxes = [];
            obj.currentFigure = [];
        end
        
        function runMouseReplicate(obj, subjectNames, funString)
            close;
            for i = 1:length(obj.subjects)
                replicatedObj = obj;
                replicatedObj.subjects = subjectNames;
                replicatedObj.blocks = repmat(replicatedObj.blocks(i), length(replicatedObj.subjects),1);
                replicatedObj.params = repmat(replicatedObj.params(i), length(replicatedObj.subjects),1);
                replicatedObj.expDate = [];
                eval(['replicatedObj.' funString ';']);
                suptitle(obj.subjects(i));
            end
            
        end
    end
    
    methods (Static)
        function [areaNames, galvoSites, laserType] = getSpecificGalvoParadigm(siteReference)
            switch lower(siteReference)
                case 'ina3'
                    areaNames = {'Left Visual';'Both Visual';'Right Visual';'Left M2';'Both M2';'Right M2';'Left Control';'Both Control';'Right Control'};
                    galvoSites = {[-2.5,-3]; [2.5,-3]; [2.5,-3]; [-1,2]; [1,2]; [1,2]; [-1.5,0]; [1.5,0]; [1.5,0]};
                    laserType = repmat([1;2;1],3,1);
            end
        end
    end
end
classdef behaviorAnalysis
    properties (Access=public)
        subjects;                %Mouse names--optional input
        expDate;                 %Recording dates to use--optional input
        blocks;
        params;
        uniqueConditions;
        xAxisLabel;
        yAxisLabel;
        axisTicks;
        gridIdx;
    end
    
    methods
        function obj = behaviorAnalysis(subjects, expDate)
            if ~exist('subjects', 'var') || isempty(subjects)
                subjects = {'PC005';'PC006';'PC010';'PC011';'PC012';'PC013';'PC015'; 'PC016'};
            end
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~iscell(subjects); subjects = {subjects}; end
            if ~iscell(expDate); expDate = {expDate}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            
            obj.subjects = subjects;
            obj.expDate = expDate;
            [obj.blocks, obj.params]  = cellfun(@(x,y) getFilesFromDates(x, y), subjects, expDate, 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
            
            retainIdx = ones(length(obj.params),1)>0;
            minNumTrials = 150;
            if any([obj.params.validTrials]<minNumTrials)
                fprintf('Warning: Removing days with less than %d trials\n', minNumTrials);
                retainIdx([obj.params.validTrials]<minNumTrials) = 0;
            end
            for i = 1:length(subjects)
                mouseIdx = strcmp({obj.blocks.subject}', subjects{i});
                [conditionSets, ~, setIdx] = unique(cellfun(@(x) num2str(x(:)'),{obj.blocks(mouseIdx).uniqueConditions}','uni',0));
                if length(conditionSets)>1; disp('Warning: Several parameter sets in date range. Using mode');
                    retainIdx(mouseIdx) = retainIdx(mouseIdx).*(setIdx == mode(setIdx));
                end
            end
            obj.blocks = obj.blocks(retainIdx);
            obj.params = obj.params(retainIdx);
            obj.subjects = unique({obj.blocks.subject}');
            [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
            obj.blocks = arrayfun(@(x) obj.blocks(subjectIdx==x), 1:length(subjects), 'uni', 0)';
            obj.params = arrayfun(@(x) obj.params(subjectIdx==x), 1:length(subjects), 'uni', 0)';
            if length(obj.subjects)~=length(subjects); disp('Warning: No valid data found for some subjects'); end
            
            for i = 1:length(subjects)
                uniqueConditions = double([diff(obj.blocks{i}(1).uniqueConditions(:,1:2),[],2) ....
                    diff(obj.blocks{i}(1).uniqueConditions(:,3:4),[],2) ...
                    obj.blocks{i}(1).uniqueConditions(:,5)]);
                obj.xAxisLabel{i,1} = 'Relative contrast';
                if size(unique(uniqueConditions, 'rows'),1) < size(obj.blocks{i}(1).uniqueConditions,1)
                    error('Unprepared to multisensory combinations of this nature');
                end
                if all([length(unique(abs(uniqueConditions(:,3)))) length(unique(abs(uniqueConditions(:,1))))])>2
                    error('Detected changes in both audAmplitude and audInitialAzimuth so cannot plot');
                end
                
                if length(unique(abs(uniqueConditions(:,3)))) > 2
                    fprintf('%s will be analyzed based on audInitialAzimuth changes', obj.subjects{i});
                    obj.yAxisLabel{i,1} = 'Auditory azimuth';
                    uniqueConditions = uniqueConditions(:,[3,2]);
                else; obj.yAxisLabel{i,1} = 'Auditory amplitude';
                    uniqueConditions = uniqueConditions(:,1:2);
                end
                obj.axisTicks{i,1} = {unique(uniqueConditions(:,1)) unique(uniqueConditions(:,2))};
                [visGridConditions, audGridConditions] = meshgrid(obj.axisTicks{i}{2}, obj.axisTicks{i}{1});
                [~, obj.gridIdx{i,1}] = ismember(uniqueConditions, [audGridConditions(:) visGridConditions(:)], 'rows');
                obj.uniqueConditions{i,1} = uniqueConditions;
            end
        end
        
        function viewBoxPlots(obj)
            figure;
            maxGrid = max(cellfun(@length, [vertcat(obj.axisTicks{:})]), [], 1);
            for i  = 1:length(obj.subjects)
                responses = vertcat(obj.blocks{i}(:).response);
                conditions = vertcat(obj.blocks{i}(:).conditions);
                fracRightTurns = nan*ones(length(obj.axisTicks{i}{1}), length(obj.axisTicks{i}{2}));
                for j = 1:length(obj.gridIdx{i})
                    fracRightTurns(obj.gridIdx{i}(j)) = mean(responses(conditions==j)==2);
                end
                obj.getAxes(obj.subjects{i}, 30, maxGrid(2)/(1.3*maxGrid(1)), [25 50], [50 100]);
                imsc(fracRightTurns, [0,1], redblue, 'k');
                daspect([1 1 1])
                title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(responses), length(obj.blocks{i})));
                set(gca, 'xTick', 1:size(fracRightTurns,2), 'xTickLabel', obj.axisTicks{i}{2}*100, 'fontsize', 14)
                set(gca, 'yTick', 1:size(fracRightTurns,1), 'yTickLabel', obj.axisTicks{i}{1}*10, 'fontsize', 14)
                xlabel(obj.xAxisLabel{i});
                ylabel(obj.yAxisLabel{i});
                box off;
            end
            currentAxisPotision = get(gca, 'position');
            colorBar = colorbar; set(colorBar,'YTick',[0,1]);
            set(gca, 'position', currentAxisPotision);
            figureSize = get(gcf, 'position');
            set(colorBar, 'position', [1-75/figureSize(3), 0.2, 30/figureSize(3), 0.6])
            colorLabel = ylabel(colorBar,'Fraction of right turns');
            set(colorLabel, 'position', [1 0.5000 0], 'FontSize', 14)
        end
        
        function viewPsychoCurves(obj)
            figure;
            maxGrid = max(cellfun(@length, [vertcat(obj.axisTicks{:})]), [], 1);
            for i  = 1:length(obj.subjects)
                responses = vertcat(obj.blocks{i}(:).response);
                conditions = vertcat(obj.blocks{i}(:).conditions);
                fracRightTurns = nan*ones(length(obj.axisTicks{i}{1}), length(obj.axisTicks{i}{2}));
                for j = 1:length(obj.gridIdx{i})
                    fracRightTurns(obj.gridIdx{i}(j)) = mean(responses(conditions==j)==2);
                end
                obj.getAxes(obj.subjects{i}, 30, maxGrid(2)/(1.3*maxGrid(1)), [25 50], [50 100]);
                imsc(fracRightTurns, [0,1], redblue, 'k');
                daspect([1 1 1])
                title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(responses), length(obj.blocks{i})));
                set(gca, 'xTick', 1:size(fracRightTurns,2), 'xTickLabel', obj.axisTicks{i}{2}*100, 'fontsize', 14)
                set(gca, 'yTick', 1:size(fracRightTurns,1), 'yTickLabel', obj.axisTicks{i}{1}*10, 'fontsize', 14)
                xlabel(obj.xAxisLabel{i});
                ylabel(obj.yAxisLabel{i});
                box off;
            end
            currentAxisPotision = get(gca, 'position');
            colorBar = colorbar; set(colorBar,'YTick',[0,1]);
            set(gca, 'position', currentAxisPotision);
            figureSize = get(gcf, 'position');
            set(colorBar, 'position', [1-75/figureSize(3), 0.2, 30/figureSize(3), 0.6])
            colorLabel = ylabel(colorBar,'Fraction of right turns');
            set(colorLabel, 'position', [1 0.5000 0], 'FontSize', 14)
        end
        
        function obj = getAxes(obj, subject, axisGap, figureRatio, botTopEdge, leftRightEdge)
            if ~exist('subject', 'var'); subject = obj.subject{1}; end
            if ~iscell(subject); subject = {subject}; end
            if ~exist('figureRatio', 'var'); figureRatio = 1; end
            if ~exist('axisGap', 'var'); axisGap = 25; end
            if ~exist('botTopEdge', 'var'); botTopEdge = 50; end
            if ~exist('leftRightEdge', 'var'); leftRightEdge = 50; end
            screenSize = get(0,'ScreenSize');
            screenRatio = round(screenSize(3)/screenSize(4));
            numOfSubjects = length(unique(obj.subjects));
            numOfRows = find(((1:5)*screenRatio.*(1:5))>=numOfSubjects,1);
            numOfCols = ceil(numOfSubjects/numOfRows);
            figureSize = min([400*numOfCols*figureRatio, 400*numOfRows], screenSize(3:4));
            
            botTopEdge = botTopEdge/figureSize(2);
            leftRightEdge = leftRightEdge/figureSize(1);
            axisGap = axisGap./figureSize;
            
            tightSubplot(numOfRows, numOfCols, find(contains(obj.subjects,subject)), axisGap, botTopEdge, leftRightEdge);
            set(gcf, 'position', [screenSize(3:4)-figureSize, figureSize])
        end
    end
end
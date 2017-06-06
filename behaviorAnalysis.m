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
                subjects = {'PC005';'PC006';'PC010';'PC011';'PC012';'PC013';'PC015'; 'PC016'}';
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
                [visGridConditions, audGridConditions] = meshgrid(obj.axisTicks{1}{1}, obj.axisTicks{1}{2});
                [~, obj.gridIdx{i,1}] = ismember(uniqueConditions, [visGridConditions(:) audGridConditions(:)], 'rows');
                
                arrayfun(@(x,y) find(all([x,y]==uniqueConditions),1), obj.gridConditions{i,1}(:,:,1), obj.gridConditions{i,1}(:,:,1))
                obj.uniqueConditions{i,1} = uniqueConditions;
                
                
%                 if length(unique(abs(obj.uniqueConditions(:,end-1:end))))
                audAmplitudeDiff = diff(vertcat(selectedBlocks.audLeftRight), [], 2);
                visContrastDiff = diff(vertcat(selectedBlocks.visLeftRight), [], 2);
                audInitialAzimuth = vertcat(selectedBlocks.audInitialAzimuth);
                visInitialAzimuth = vertcat(selectedBlocks.visInitialAzimuth);
                
                uniqueAudAmplitudeDiff = unique(audAmplitudeDiff);
                uniqueAudInitialAzimuth = unique(audInitialAzimuth);
                uniqueVisContrastDiff = unique(visContrastDiff);
                uniqueVisInitialAzimuth = unique(visInitialAzimuth);
                
                if length(uniqueAudInitialAzimuth)>2 && length(uniqueAudAmplitudeDiff)>3
                    fprintf('Warning: skipping %s because too many audio parameters', obj.subjects{i});
                elseif length(uniqueAudInitialAzimuth)>2
                    yAxisLabel = 'Audio Azimuth'; uniqueAudConditions = uniqueVisInitialAzimuth;
                    audConditions(:,1) = uniqueConditions(:,5);
                elseif length(uniqueAudAmplitudeDiff)>3
                    yAxisLabel = 'Audio Amplitude'; uniqueAudConditions = uniqueAudAmplitudeDiff;
                    audConditions(:,1) = diff(uniqueConditions(:,1:2),[],2);
                end
                
                if length(uniqueVisInitialAzimuth)>2 && length(uniqueVisContrastDiff)>3
                    fprintf('Warning: skipping %s because too many visual parameters', obj.subjects{i});
                elseif length(uniqueVisInitialAzimuth)>2
                    xAxisLabel = 'Visual Contrast'; uniqueVisConditions = uniqueVisInitialAzimuth;
                    visConditions(:,1) = uniqueConditions(:,5);
                elseif length(uniqueVisContrastDiff)>3
                    xAxisLabel = 'Visual Contrast'; uniqueVisConditions = uniqueVisContrastDiff;
                    visConditions(:,1) = diff(uniqueConditions(:,1:2),[],2);
                end
            end
        end
        
        function viewBoxPlots(obj)
            for i  = 1:length(obj.subjects)
                selectedBlocks = obj.blocks(strcmp({obj.blocks.subject}', obj.subjects{i}));
                uniqueConditions = selectedBlocks(1).uniqueConditions;
                audVis = unique(diff(vertcat(obj.blocks.audLeftRight)));
                
                
                
            end
            rTyp = unique(cat(1,aBlk.cRes));
            if length(rTyp) > 2; error('Too many response types'); end
            [rMat, aLab] = getRespMatrix(aBlk);
            
            axes(curA);
            if min(size(rMat(:,:,1,1))) == 1; imagesc(mean(rMat(:,:,:,1),3)); ...
                    colormap(redblue); caxis([0 1]);
            else; try imsc(mean(rMat(:,:,:,1),3), [0,1], redblue, 'k');
                catch; imagesc(mean(rMat(:,:,:,1),3)); end
            end
            warning('off', 'all'); axis equal tight; box off; warning('on', 'all');
            %%
            set(gca, 'xTick', 1:size(rMat,2), 'xTickLabel', aLab{2}*100, 'fontsize', 14)
            set(gca, 'yTick', 1:size(rMat,1), 'yTickLabel', aLab{1}, 'fontsize', 14)
            cBar = colorbar; set(cBar,'YTick',[0,1])
            title(sprintf('%s: %d trials, %d sessions', aBlk(1).mNam, length(cat(1,aBlk.cRes)), length(aBlk)));
            % xlabel('Contrast Difference (%)');
            % ylabel('Audio Amplitude Difference');
            % boxM(isnan(boxM)) = -0.1;
        end
        
        function obj = getAxes(obj, subject)
            if ~exist('subject', 'var'); subject = obj.subject{1}; end
            if ~iscell(subject); subject = {subject}; end
            screenRatio = get(0,'ScreenSize');
            screenRatio = round(screenRatio(3)/screenRatio(4));
            numOfSubjects = length(unique(obj.subject));
            numOfRows = find(((1:5)*screenRatio.*(1:5))>=numOfSubjects,1);
            numOfCols = ceil(numOfSubjects/numOfRows);
            tightSubplot(numOfRows,numOfCols, find(contains(obj.subject,subject)), 0.02);
        end
    end
end
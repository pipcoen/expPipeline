classdef GLMmulti < matlab.mixin.Copyable
    properties (Access=public)
        modelString;
        prmLabels;
        prmFits;
        prmBounds;
        prmInit;
        inputFun;
        modelFun;
        blockData;
        pHat;
        evalPoints;
        initGuess;
    end
    
    methods
        function obj = GLMmulti(inputblockData)
            %% Input blockData must be a struct with fields: conditions and responseMade
            obj.blockData = inputblockData;
            tab = tabulate(obj.blockData.responseMade)/100;
            obj.initGuess = sum(tab(:,3).*log2(tab(:,3)));
        end
        
        function obj = setModel(obj,modelString,initParams)
            % Function which contains the model definitions. Each definition requires a set of parameter labels, estimation bounds, and a (usually)
            % anonymous function. The function (modelFun)define the linear model in 2-class multinomial regression for 2-choice behaviour. The
            % inputFun function provides input blockData to modelFun.
            %
            % modelFun usually takes two inputs P and IN. P is a vector of param values. IN is a matrix of model inputs derived from inputFun
            %
            % inputFun usually takes in the blockData struct and outputs a matrix of variables to be used in the mode;  
                
            obj.modelString = modelString;            
            obj.blockData.audValues = unique(obj.blockData.audDiff);
            obj.blockData.visValues = unique(obj.blockData.visDiff);
            
            uniA = obj.blockData.audValues;
            uniV = obj.blockData.visValues;
            comb = unique([obj.blockData.visDiff obj.blockData.audDiff], 'rows');
            numA = length(uniA);
            numV = length(uniV);
            numC = length(comb);
            
            
            maxContrast = max(abs(uniV));
            audRepMat = repmat(uniA,1,200)';
            obj.evalPoints = [repmat(linspace(-maxContrast,maxContrast,200)', numA,1), audRepMat(:)];
            [aGrid,vGrid] = meshgrid(uniA,uniV);
            
            %%
            indiSum = @(pRange,diffVals) ...
                sum(cell2mat(arrayfun(@(x,y,z) x.*(all(y{1}==z{1},2)), pRange,repmat({diffVals},1,length(unique(diffVals, 'rows'))), num2cell(unique(diffVals, 'rows'),2)', 'uni', 0)),2);
            
            sqrtDiff = @(diffVals) ...
                sqrt(abs(diffVals)).*sign(diffVals);
            
            switch(modelString)
                case 'SimpleLogistic'
                    obj.prmLabels = {'bias','visScale','audScale'};
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)([b.visDiff, b.audDiff]);
                    obj.modelFun = @(P,in) (P(1)+ P(2)*in(:,1) + P(3)*in(:,2));
                case 'SqrtLogistic'
                    obj.prmLabels = [{'bias','visScale'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)'];
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)([b.visDiff, b.audDiff]);
                    obj.modelFun = @(P,in) (P(1)+ P(2)*sqrtDiff(in(:,1))+indiSum(P(3:(2+numA)),in(:,2)));
                case 'SqrtLogisticSplit'
                    obj.prmLabels = [{'bias','visScaleL','visScaleR'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)'];
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)([b.visDiff, b.audDiff]);
                    obj.modelFun = @(P,in) (P(1)+ P(2)*sqrtDiff(in(:,1).*(in(:,1)>0))+P(3)*sqrtDiff(in(:,1).*(in(:,1)<0))+indiSum(P(4:(3+numA)),in(:,2)));
                case 'SqrtLogisticSplitDelta'
                    obj.prmLabels = [{'bias','visScaleL','visScaleR'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)'];
                    obj.prmLabels = cellfun(@(x) [x 'Delta'], obj.prmLabels, 'uni', 0);
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)({b.visDiff; b.audDiff; b.prmInit});
                    obj.modelFun = @(P,in) (P(1)+in{3}(1) + (P(2)+in{3}(2))*sqrtDiff(in{1}.*(in{1}>0)) + (P(3)++in{3}(3))*sqrtDiff(in{1}.*(in{1}<0)) + ...
                        indiSum(P(4:(3+numA))+in{3}(4:(3+numA)),in{2}));
                case 'Simp-emp'
                    allValues = [cellfun(@(x) [num2str(x) 'Vis'], num2cell(uniV), 'uni', 0); cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)];
                    obj.prmLabels = ['bias'; allValues(:)]';
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)([b.visDiff, b.audDiff]);
                    obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numV)),in(:,1)) + indiSum(P((2+numV):(1+numV+numA)),in(:,2)));
                    obj.evalPoints = [vGrid(:) aGrid(:)];
                case 'Simp-aud'
                    allValues = cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0);
                    obj.prmLabels = ['bias'; allValues(:)]';
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)([b.visDiff, b.audDiff]);
                    obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numA)),in(:,2)));
                    obj.evalPoints = [vGrid(:) aGrid(:)];
                case 'Simp-vis'
                    allValues = cellfun(@(x) [num2str(x) 'Vis'], num2cell(uniV), 'uni', 0);
                    obj.prmLabels = ['bias'; allValues(:)]';
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)([b.visDiff, b.audDiff]);
                    obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numV)),in(:,1)));
                    obj.evalPoints = [vGrid(:) aGrid(:)];
                case 'Full-emp'
                    allValues = [cellfun(@(x) [sprintf('%0.1f', x), 'VisAud'], num2cell(comb,2), 'uni', 0)];
                    obj.prmLabels = ['bias'; allValues]';
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.inputFun = @(b)([b.visDiff, b.audDiff]);
                    obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numC)),in));
                    obj.evalPoints = comb;
                otherwise
                    error('Model does not exist');
            end
            
            if isempty(obj.prmInit)
                obj.prmInit = zeros(1,length(obj.prmLabels));
            end
        end
        
        function obj = fit(obj)
            %Non crossvalidated fitting
            if isempty(obj.modelFun); error('Set model first'); end
            options = optimset('algorithm','interior-point','MaxFunEvals',100000,'MaxIter',10000);
            fittingObjective = @(b) (obj.calculateLogLik(b, obj.inputFun(obj.blockData), obj.blockData.responseMade));
            [obj.prmFits,~,exitflag] = fmincon(fittingObjective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
            if ~any(exitflag == [1,2])
                obj.prmFits = nan(1,length(obj.prmLabels));
            end
        end
        
        function [obj,logLik] = fitCV(obj,nFolds)
            %Crossvalidated fitting
            if isempty(obj.modelFun); error('Set model first'); end
            if ~exist('nFolds', 'var') || isempty(nFolds); nFolds = length(obj.blockData.responseMade); end
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',2000);
            cvObj = cvpartition(obj.blockData.responseMade,'KFold',nFolds);
            
            obj.prmFits = nan(cvObj.NumTestSets,length(obj.prmLabels));
            obj.pHat = [];
            allInputs = obj.inputFun(obj.blockData);
            logLik = nan(cvObj.NumTestSets,1);
            for i = 1:cvObj.NumTestSets
                disp(['Model: ' obj.modelString '. Fold: ' num2str(i) '/' num2str(cvObj.NumTestSets)]);
                
                fittingObjective = @(b) (obj.calculateLogLik(b,  allInputs(cvObj.training(i),:), obj.blockData.responseMade(cvObj.training(i))));
                [obj.prmFits(i,:),~,exitflag] = fmincon(fittingObjective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
                if ~any(exitflag == [1,2]); obj.prmFits(i,:) = nan(1,length(obj.prmLabels)); end
                
                pHatTrained = obj.calculatepHat(obj.prmFits(i,:), allInputs(cvObj.test(i),:));
                
                testIdx = sub2ind(size(pHatTrained), 1:size(pHatTrained, 1), obj.blockData.responseMade(cvObj.test(i))');
                obj.pHat(cvObj.test(i)) = pHatTrained(testIdx);
                logLik(i)=mean(-log2(obj.pHat(cvObj.test(i))));
            end
        end
        
        function h = plotblockData(obj)
            %%
            numTrials = prc.makeGrid(obj.blockData, obj.blockData.responseMade~=0, @length, 1);
            numRightTurns = prc.makeGrid(obj.blockData, obj.blockData.responseMade==2, @sum, 1);
            
            audValues = obj.blockData.audValues;
            colorChoices = plt.selectRedBlueColors(audValues);
            
            [prob,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), numRightTurns, numTrials, 'uni', 0);
            prob = cell2mat(cellfun(@(x) permute(x, [3,1,2]), prob, 'uni', 0));
            lowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
            highBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));
            
            for audVal = audValues(:)'
                idx = find(obj.blockData.grids.audValues==audVal & numTrials>0);
                err = [prob(idx)-lowBound(idx), highBound(idx) - prob(idx)];
                errorbar(obj.blockData.grids.visValues(idx),prob(idx),err(:,1),err(:,2),'.','MarkerSize',20, 'Color', colorChoices(audValues==audVal,:));
                hold on;
            end
            xlabel('Contrast');
            ylabel('P( choice | contrast)');
            set(gca,'box','off');
            h=gca;
            set(gcf,'color','w');
        end
        
        function figureHand = plotFit(obj)
            if size(obj.prmFits,1)~=1; error('Model not fitted (non-crossvalidated) yet'); end
            hold on;
            colorChoices = plt.selectRedBlueColors(obj.blockData.audValues);
            otherInputs = obj.inputFun(obj.blockData);
            otherInputs(:,1:2)=[];
            if isempty(otherInputs); inputs = obj.evalPoints;
            else, inputs = [obj.evalPoints, zeros(length(obj.evalPoints),size(otherInputs,2))];
            end
            pHatCalculated = obj.calculatepHat(obj.prmFits,inputs);
            for audVal = obj.blockData.audValues(:)'
                plotIdx = obj.evalPoints(:,2)==audVal;
                plot(gca, obj.evalPoints(plotIdx,1), pHatCalculated(plotIdx,2), 'Color', colorChoices(obj.blockData.audValues==audVal,:), 'linewidth', 2);
            end
            maxContrast = max(abs(obj.evalPoints(:,1)));
            xlim([-maxContrast maxContrast])
            set(gca, 'xTick', (-maxContrast):(maxContrast/4):maxContrast, 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)*100));
            title(obj.modelString);
            hold off;
            figureHand = gca;
        end
        
        
        function plotPredVsActual(obj, varargin)
            
            plotParams.Color = 'r';
            plotParams.ax = gca;
            if ~isempty(varargin)
                plotParams = mergeStructs(varargin{1},plotParams);
            end
            
            switch(obj.ContrastDimensions)
                case 1
                    contrast1D = diff(obj.blockData.stimulus, [], 2);
                    uniqueC1D = unique(contrast1D);
                    nC = length(uniqueC1D);
                    prop=zeros(nC,3);
                    prop_ci=zeros(nC,3,2);
                    for c = 1:length(uniqueC1D)
                        D = obj.getrow(obj.blockData,contrast1D == uniqueC1D(c));
                        respSum = sum([D.responseMade==1 D.responseMade==2 D.responseMade==3],1);
                        p = respSum/length(D.responseMade);
                        
                        [p,pci]=binofit(respSum,length(D.responseMade),0.05);
                        
                        prop_ci(c,:,:) = pci;
                        
                        prop(c,:) = p;
                    end
                    
                    if max(obj.blockData.responseMade) == 2 %for 2AFC tasks
                        rMax = 2;
                    else
                        rMax = 3;
                    end
                    
                    evalCon = unique(obj.blockData.stimulus,'rows');
                    evalC1d = evalCon(:,2) - evalCon(:,1);
                    [~,sortIdx]=sort(evalC1d);
                    evalCon = evalCon(sortIdx,:);
                    pHat = obj.calculatepHat(obj.prmFits,evalCon);
                    
                    rSymbols = {'o', '.', 'x'};
                    
                    for c = 1:length(uniqueC1D)
                        for r = 1:rMax
                            plot(plotParams.ax, pHat(c,r), prop(c,r), rSymbols{r}, 'Color', plotParams.Color)
                            hold on;
                        end
                        for r = 1:rMax
                            plot(plotParams.ax, pHat(c,r)*ones(1,2), squeeze(prop_ci(c,r,:)), 'Color', plotParams.Color)
                        end
                    end
                    plot(plotParams.ax, [0 1], [0 1], 'k--');
                    
                    ylabel('actual probability');
                    xlabel('predicted probability');
                    legend({'left resp', 'right resp', 'nogo'}, 'Location', 'Best');
                    
                    axis square
                    box off
                    
                case 2
                    fprintf(1, 'plotPredVsActual not yet implemented for 2D task\n')
            end
        end
        
        function h = plotParams(obj)
            if size(obj.prmFits,1)~=1; return; end
            bar(obj.prmFits);
            set(gca,'XTickLabel',obj.prmLabels,'XTick',1:numel(obj.prmLabels));
            title(obj.modelString);
            h=gca;
        end
        
        function pHat = calculatepHat(obj,testParams,inputs)
            if isempty(obj.modelFun); error('Set model first'); end
            logOdds = obj.modelFun(testParams,inputs);
            probLeft = exp(logOdds)./(1+exp(logOdds));
            pHat = [probLeft 1-probLeft zeros(length(probLeft),1)];
        end
        
        function logLik = calculateLogLik(obj,testParams,inputs,responses)
            pHatCalculated = obj.calculatepHat(testParams, inputs);
            logLik = -sum(log2( pHatCalculated(sub2ind(size(pHatCalculated), (1:length(responses))', responses)) ));
        end
    end
end
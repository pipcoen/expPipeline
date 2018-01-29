classdef GLMmulti
    properties (Access=public)
        modelString;
        prmLabels;
        prmFits;
        prmBounds;
        prmInit;
        Zinput;
        ZL;
        ZR;
        data;
        p_hat;
        evalPoints;
    end
    
    properties (Access=public)
        guess_bpt;
        lapseFlag=0;
    end
    
    methods
        function obj = GLMmulti(inputData)
            %% Input data must be a struct with fields: conditions and response
            obj.data = inputData;
            
            tab = tabulate(obj.data.responseMade);
            tab = tab(:,3)/100;
            obj.guess_bpt=sum(tab.*log2(tab));
        end
        
        function obj = setModel(obj,modelString)
            % Function which contains the model definitions. Each definition requires a set of prm Labels, estimation
            % bounds, and some anonymous functions. The functions ZL and ZR define the two linear models in 3-class multinomial
            % regression for 3-choice behaviour. For 2-choice behaviour, only ZL needs to be defined for logistic regression. The
            % Zinput function provides input data to ZL and ZR. Often this will be simply the trial contrast but it can also be
            % any other feature of interest.
            %
            % ZL and ZR always take two inputs P and IN. P is a vector of prm values. IN is a matrix of model inputs as derived
            % from the Zinput function
            %
            % Zinput always takes in the data struct and outputs a matrix of variables to be used in the model
            
            obj.modelString = modelString;
            obj.prmFits = [];
            obj.prmInit = [];
            
            obj.data.audValues = unique(obj.data.audDiff);
            obj.data.visValues = unique(obj.data.visDiff);
            maxContrast = max(abs(obj.data.visValues));
            audRepMat = repmat(obj.data.audValues,1,200)';
            obj.evalPoints = [repmat(linspace(-maxContrast,maxContrast,200)', length(obj.data.audValues),1), audRepMat(:)];
            switch(modelString)
                case 'SimpleLogistic'
                    obj.prmLabels = {'bias','visScale','audScale'};
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.Zinput = @(b)([b.visDiff, b.audDiff]);
                    obj.ZL = @(P,in) (P(1)+ P(2)*in(:,1) + P(3)*in(:,2));
                    obj.ZR = [];
                case 'Simp-emp'
                    allValues = [cellfun(@(x) [num2str(x) 'Vis'], num2cell(obj.data.visValues), 'uni', 0); cellfun(@(x) [num2str(x) 'Aud'], num2cell(obj.data.audValues), 'uni', 0)];              
                    obj.prmLabels = ['bias'; allValues]';
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.Zinput = @(b)([b.visDiff, b.audDiff]);
                    obj.ZL = @(P,in)(P(1) +  ...
                        sum(cell2mat(arrayfun(@(x,y,z) x.*(y{1}==z), P(2:(length(unique(in(:,1)))+1)),repmat({in(:,1)},1,length(unique(in(:,1)))), unique(in(:,1))', 'uni', 0)),2) + ...
                        sum(cell2mat(arrayfun(@(x,y,z) x.*(y{1}==z), P((length(unique(in(:,1)))+2):(length(unique(in(:,1)))+length(unique(in(:,2)))+1)),repmat({in(:,2)},1,length(unique(in(:,2)))), unique(in(:,2))', 'uni', 0)),2));
                    obj.ZR = [];
                    [aGrid,vGrid] = meshgrid(obj.data.audValues,obj.data.visValues);
                    obj.evalPoints = [vGrid(:) aGrid(:)];                    
                otherwise
                    error('Model does not exist');
            end
            
            if isempty(obj.prmInit)
                obj.prmInit = zeros(1,length(obj.prmLabels));
            end
        end
        
        function obj = addLapse(obj)
            if isempty(obj.ZL); error('Set model first'); end
            obj.prmLabels = [obj.prmLabels 'LapseRate'];
            obj.prmBounds = [obj.prmBounds, [0;1]];
            obj.prmInit = [obj.prmInit 0];
        end
        
        function obj = fit(obj)
            %Non crossvalidated fitting
            if isempty(obj.ZL); error('Set model first'); end
            options = optimset('algorithm','interior-point','MaxFunEvals',100000,'MaxIter',10000);
            fittingObjective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), obj.data.responseMade));
            [obj.prmFits,~,exitflag] = fmincon(fittingObjective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
            if ~any(exitflag == [1,2])
                obj.prmFits = nan(1,length(obj.prmLabels));
            end
        end
        
        function r2(obj)
            if isempty(obj.prmFits)
                error('Fit first');
            end
            
            %non-cv log likelihood for fitted model
            ll_model = obj.calculateLogLik(obj.prmFits,obj.data.stimulus,obj.data.responseMade)/length(obj.data.responseMade);
            
            %another thing to do is just take the mode of the pred
            %probabilities and see if they match the data
            pHats = obj.calculatepHat(obj.prmFits,obj.data.stimulus);
            classify = nan(length(obj.data.responseMade),1);
            for t = 1:length(obj.data.responseMade)
                [~,rhat] = max(pHats(t,:));
                
                if rhat == obj.data.responseMade(t)
                    classify(t)=1;
                else
                    classify(t)=0;
                end
            end
                      
            r2 = 1 + ll_model/obj.guess_bpt;
            disp(['McFadden Pseudo-R^2: ' num2str(r2)]);
            disp(['Proportion correctly classified: ' num2str(mean(classify))]);
            
            
        end
        
        function [obj,varargout] = fitCV(obj,varargin)
            %Crossvalidated fitting
            
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',2000);
            
            if isempty(varargin)
                C = cvpartition(length(obj.data.responseMade),'LeaveOut');
            else
                C = cvpartition(obj.data.responseMade,'KFold',varargin{1});
            end
            
            obj.prmFits = nan(C.NumTestSets,length(obj.prmLabels));
            obj.p_hat = [];
            for f=1:C.NumTestSets
                disp(['Model: ' obj.modelString '. Fold: ' num2str(f) '/' num2str(C.NumTestSets)]);
                trainIdx = find(C.training(f)==1);
                testIdx = find(C.test(f)==1);
                
                inputs = obj.Zinput(obj.data);
                trainInputs = inputs(trainIdx,:);
                testInputs = inputs(testIdx,:);
                
                trainResponses = obj.data.responseMade(trainIdx);
                testResponse = obj.data.responseMade(testIdx);
                
                objective = @(b) ( obj.calculateLogLik(b, trainInputs, trainResponses) );

                [obj.prmFits(f,:),~,exitflag] = fmincon(objective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
                if ~any(exitflag == [1,2])
                    obj.prmFits(f,:) = nan(1,length(obj.prmLabels));
                end
                %                 end
                
                
                pHat = obj.calculatepHat(obj.prmFits(f,:), testInputs);
                
                for i = 1:length(testResponse)
                    obj.p_hat(testIdx(i),1) = pHat(i,testResponse(i));
                end
                
                LL(f)=mean(-log2(obj.p_hat(testIdx)));
            end
            
            varargout = {LL};
        end
        
        function h = plotData(obj)
            %%
            numTrials = prc.makeGrid(obj.data, obj.data.responseMade~=0, @length, 1);
            numRightTurns = prc.makeGrid(obj.data, obj.data.responseMade==2, @sum, 1);
            
            audValues = obj.data.audValues;
            colorChoices = plt.selectRedBlueColors(audValues);
            
            [prob,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), numRightTurns, numTrials, 'uni', 0);
            prob = cell2mat(cellfun(@(x) permute(x, [3,1,2]), prob, 'uni', 0));
            lowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
            highBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));
            
            for audVal = audValues(:)'
                idx = find(obj.data.grids.audValues==audVal & numTrials>0);
                err = [prob(idx)-lowBound(idx), highBound(idx) - prob(idx)];
                errorbar(obj.data.grids.visValues(idx),prob(idx),err(:,1),err(:,2),'.','MarkerSize',20, 'Color', colorChoices(audValues==audVal,:));
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
            colorChoices = plt.selectRedBlueColors(obj.data.audValues);
            otherInputs = obj.Zinput(obj.data);
            otherInputs(:,1:2)=[];
            if isempty(otherInputs); inputs = obj.evalPoints;
            else, inputs = [obj.evalPoints, zeros(length(obj.evalPoints),size(otherInputs,2))];
            end
            pHat = obj.calculatepHat(obj.prmFits,inputs);            
            for audVal = obj.data.audValues(:)'
                plotIdx = obj.evalPoints(:,2)==audVal;
                plot(gca, obj.evalPoints(plotIdx,1), pHat(plotIdx,2), 'Color', colorChoices(obj.data.audValues==audVal,:), 'linewidth', 2);
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
                    contrast1D = diff(obj.data.stimulus, [], 2);
                    uniqueC1D = unique(contrast1D);
                    nC = length(uniqueC1D);
                    prop=zeros(nC,3);
                    prop_ci=zeros(nC,3,2);
                    for c = 1:length(uniqueC1D)
                        D = obj.getrow(obj.data,contrast1D == uniqueC1D(c));
                        respSum = sum([D.responseMade==1 D.responseMade==2 D.responseMade==3],1);
                        p = respSum/length(D.responseMade);
                        
                        [p,pci]=binofit(respSum,length(D.responseMade),0.05);
                        
                        prop_ci(c,:,:) = pci;
                        
                        prop(c,:) = p;
                    end
                    
                    if max(obj.data.responseMade) == 2 %for 2AFC tasks
                        rMax = 2;
                    else
                        rMax = 3;
                    end
                    
                    evalCon = unique(obj.data.stimulus,'rows');
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
            if size(obj.prmFits,1)==1
                bar(obj.prmFits);
                set(gca,'XTickLabel',obj.prmLabels,'XTick',1:numel(obj.prmLabels));
                title(obj.modelString);
                h=gca;
            end
        end
        
        function pHat = calculatepHat(obj,testParams,inputs)
            if isempty(obj.ZL); error('Set model first'); end

            switch(obj.lapseFlag)
                case 0
                    lapse = 0;
                case 1
                    lapse = testParams(end);
                case 2
                    lapse = testParams(end);
            end
            
            if isempty(obj.ZR) %if a AFC task then no ZR is defined, only pL vs pR
                zl = obj.ZL(testParams,inputs);
                pL = lapse + (1-2*lapse)*exp(zl)./(1+exp(zl));
                pR = 1 - pL;
                N = length(pL);
                pHat = [pL pR zeros(N,1)];
            else
                zl = obj.ZL(testParams,inputs);
                zr = obj.ZR(testParams,inputs);
                pL = (1-lapse)*exp(zl)./(1+exp(zl)+exp(zr));
                pR = (1-lapse)*exp(zr)./(1+exp(zl)+exp(zr));
                pNG = 1 - (pL + pR);
                
                pHat = [pL pR pNG];
            end
        end
        
        function cov=paramCov(obj,varargin)
            if isempty(obj.prmFits)
                error('Fit model first!');
            end
            
            c = obj.data.stimulus;
            H = obj.hessian(obj.prmFits,c);
            
            %Calculate fisher info
            F = -sum(H,3);
            cov = inv(F);
            
            figure;
            ax1=subplot(1,2,1);
            imagesc(cov); title('Covariance'); axis square;
            ax2=subplot(1,2,2);
            imagesc(corrcov(cov)); caxis([-1 1]); title('Correlation'); axis square;
            set(get(gcf,'children'),'XTickLabel',{'B_L','B_R','S_L','S_R'},'xtick',1:4,'yTickLabel',{'B_L','B_R','S_L','S_R'},'ytick',1:4);
            
            cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
            colormap(ax2,cmap); colorbar;
        end
        
    end
    
    
    
    methods (Access= {?GLM})
        function logLik = calculateLogLik(obj,testParams,inputs,responses)
            pHat = obj.calculatepHat(testParams, inputs);
            logLik = -sum(log2( pHat(sub2ind(size(pHat), (1:length(responses))', responses)) ));
        end
        
        function H = hessian(obj,params,contrasts)
            
            %Use untransformed contrasts in ZL and ZR eqns
            zl=obj.ZL(params,contrasts);
            zr=obj.ZR(params,contrasts);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            
            n = obj.prmFits(end-1);
            c50 = obj.prmFits(end);
            cfn = @(c)(c.^n)./(c.^n + c50.^n);
            CL = cfn(contrasts(:,1));
            CR = cfn(contrasts(:,2));
            
            H = zeros(4,4,length(contrasts));
            H(1,1,:) = pL.*(pL-1);
            H(2,2,:) = pR.*(pR-1);
            H(3,3,:) = (CL.^2).*pL.*(pL-1);
            H(4,4,:) = (CR.^2).*pR.*(pR-1);
            H(2,1,:) = pL.*pR;
            H(1,2,:) = pL.*pR;
            H(3,1,:) = CL.*pL.*(pL-1);
            H(1,3,:) = CL.*pL.*(pL-1);
            H(4,1,:) = CR.*pL.*pR;
            H(1,4,:) = CR.*pL.*pR;
            H(3,2,:) = CL.*pL.*pR;
            H(2,3,:) = CL.*pL.*pR;
            H(4,2,:) = CR.*pR.*(pR-1);
            H(2,4,:) = CR.*pR.*(pR-1);
            H(4,3,:) = CL.*CR.*pL.*pR;
            H(3,4,:) = CL.*CR.*pL.*pR;
        end
    end
end
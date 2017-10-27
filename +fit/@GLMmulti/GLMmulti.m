classdef GLMmulti
    properties (Access=public)
        modelString;
        prmLabels;
        prmFits;
        prmBounds;
        prmInit;
        prmSenses;
        Zinput;
        ZL;
        ZR;
        data;
        p_hat;
    end
    
    properties (Access=public)
        guess_bpt;
        lapseFlag=0;
    end
    
    methods
        function obj = GLMmulti(inputData)
            %% Input data must be a struct with fields: conditions and response
            obj.data = inputData;
            
            tab = tabulate(obj.data.response);
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
            
            switch(modelString)
                case 'Offset' %Model guesses based on the proportion of responses in the data
                    %used as a baseline to compare other models
                    obj.prmLabels = {'Offset_L','Offset_R'};
                    obj.prmBounds = [-inf -inf; +inf +inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1)*ones(length(in(:,1)),1));
                    obj.ZR = @(P,in)(P(2)*ones(length(in(:,1)),1));
                case 'C-subset-visual'
                    obj.prmLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R'};
                    obj.prmBounds = [-inf -inf -inf -inf;
                        +inf +inf +inf +inf];
                    obj.Zinput = @(b)([b.visContrast.*b.visInitialAzimuth<0 b.visContrast.*b.visInitialAzimuth>0]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1));
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2));
                case 'C-subset-multiAud1'
                    obj.prmLabels = {'Bias','visScale','audScale'};
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.Zinput = @(b)(double([b.visDiff sign(b.audDiff)]));
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1) + P(3).*in(:,2));
                    obj.ZR = [];
                case 'C-subset-multiAud2'
                    obj.prmLabels = {'visBias','audBias','visScale','audScale', 'bias'};
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.Zinput = @(b)(double([b.visDiff b.audDiff]));
                    obj.ZL = @(P,in)(P(1).*double(in(:,1)==0) + P(2).*double(in(:,2)==0) + P(3).*in(:,1) + P(4).*in(:,2) + P(5));
                    obj.ZR = [];
                case 'C-subset-multiAud3'
                    obj.prmLabels = {'visBias','audBias','visLScale','visRScale','audLScale','audRScale', 'bias'};
                    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
                    obj.Zinput = @(b)([b.visContrastLeftRight, b.audAzimuthLeftRight]);
                    obj.ZL = @(P,in)(P(1).*double(all(in(:,1:2),2)==0) + P(2).*double(all(in(:,3:4),2)==0) + P(3).*in(:,1) + P(4).*in(:,2) + P(5).*in(:,3) + P(6).*in(:,4) + P(7));
                    obj.ZR = [];
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

            fittingObjective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), obj.data.response));
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
            ll_model = obj.calculateLogLik(obj.prmFits,obj.data.stimulus,obj.data.response)/length(obj.data.response);
            
            %another thing to do is just take the mode of the pred
            %probabilities and see if they match the data
            pHats = obj.calculatepHat(obj.prmFits,obj.data.stimulus);
            classify = nan(length(obj.data.response),1);
            for t = 1:length(obj.data.response)
                [~,rhat] = max(pHats(t,:));
                
                if rhat == obj.data.response(t)
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
                C = cvpartition(length(obj.data.response),'LeaveOut');
            else
                C = cvpartition(obj.data.response,'KFold',varargin{1});
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
                
                trainResponses = obj.data.response(trainIdx);
                testResponse = obj.data.response(testIdx);
                
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
            numTrials = prc.makeGrid(obj.data, obj.data.response~=0, @length, 1);
            numRightTurns = prc.makeGrid(obj.data, obj.data.response==2, @sum, 1);
            
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
            h=obj.plotData();
            if size(obj.prmFits,1)~=1; error('Model not fitted (non-crossvalidated) yet'); end
            hold on;
            audValues = obj.data.audValues;
            colorChoices = plt.selectRedBlueColors(audValues);
            
            for audVal = audValues(:)'
                maxContrast = max(obj.data.visContrastLeftRight(:));
                audInputs = unique(obj.data.audAzimuthLeftRight(obj.data.audDiff==audVal,:), 'rows');
                evalPoints = [[linspace(maxContrast,0,100)'; zeros(100,1)], [zeros(100,1); linspace(0,maxContrast,100)'], repmat(audInputs, 200,1)];
                otherInputs = obj.Zinput(obj.data);
                otherInputs(:,1:4)=[];
                if isempty(otherInputs); inputs = evalPoints;
                else, inputs = [evalPoints, zeros(length(evalPoints),size(otherInputs,2))];
                end
                
                pHat = obj.calculatepHat(obj.prmFits,inputs);
                plot(h, evalPoints(:,2)-evalPoints(:,1), pHat(:,2), 'Color', colorChoices(audValues==audVal,:));
            end
            title(obj.modelString);
            hold off;
            figureHand = gca;
        end
        
        function plotPedestal(obj)
            if isempty(obj.prmFits)
                error('Need to fit model first');
            end
            
            cVals = unique(obj.data.stimulus(:));
            prop=nan(length(cVals),length(cVals),3);
            
            for cl = 1:length(cVals)
                for cr = 1:length(cVals)
                    E = obj.getrow(obj.data,obj.data.stimulus(:,1) == cVals(cl) & obj.data.stimulus(:,2) == cVals(cr));
                    for i=1:3
                        prop(cl,cr,i) = sum(E.response==i)/length(E.response);
                    end
                    pd.propChooseLeft(cl,cr) = prop(cl,cr,1);
                    pd.propChooseRight(cl,cr) = prop(cl,cr,2);
                    pd.propChooseNogo(cl,cr) = prop(cl,cr,3);
                end
            end
            
            cTestVals = min(cVals):0.01:2;
            selectContrast{1} = [reshape(reshape(repmat(cVals, 1, length(cTestVals)), length(cVals), length(cTestVals))', 1, length(cVals)*length(cTestVals)); ...
                repmat(cTestVals, 1, length(cVals))];
            selectContrast{2} = selectContrast{1}([2 1], :);
            
            predictionsSelect{1} = obj.calculatepHat(obj.prmFits, selectContrast{1}')';
            predictionsSelect{2} = obj.calculatepHat(obj.prmFits, selectContrast{2}')';
            f3 = figure; %set(f3, 'Position', [  -1896         507        1058         405]);
            
            colors = [        0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250];
            
            for ped = 1:length(cVals)-1
                
                subplot(1, length(cVals)-1, ped)
                
                for r = 1:3
                    plot(-(cVals(ped:end)-cVals(ped)), prop(ped:end,ped,r), 'o','Color',colors(r,:));
                    
                    hold on;
                    plot((cVals(ped:end)-cVals(ped)), prop(ped,ped:end,r), 'o','Color',colors(r,:));
                    
                    
                    
                    plot(-(cTestVals(cTestVals>=cVals(ped))-cVals(ped)), predictionsSelect{2}(r,selectContrast{1}(1,:)==cVals(ped)&selectContrast{1}(2,:)>=cVals(ped)),'Color',colors(r,:));
                    
                    plot((cTestVals(cTestVals>=cVals(ped))-cVals(ped)), predictionsSelect{1}(r,selectContrast{2}(2,:)==cVals(ped)&selectContrast{2}(1,:)>=cVals(ped)),'Color',colors(r,:));
                    
                end
                
                title(['pedestal = ' num2str(cVals(ped))]);
                xlabel('delta C');
                ylim([0 1]);
                xlim([-1 1]);
                set(gca,'box','off');
                %                 makepretty
            end
            set(gcf,'color','w');
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
                        respSum = sum([D.response==1 D.response==2 D.response==3],1);
                        p = respSum/length(D.response);
                        
                        [p,pci]=binofit(respSum,length(D.response),0.05);
                        
                        prop_ci(c,:,:) = pci;
                        
                        prop(c,:) = p;
                    end
                    
                    if max(obj.data.response) == 2 %for 2AFC tasks
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
        
        function obj = addData(obj,type)
            %Adds extra data depending on what is required
            
            if strcmp(obj.expRef,'none')
                error('not coded for custom struct inputs');
            end
            
            switch(type)
                case 'lick'
                    block = dat.loadBlock(obj.expRef);
                    trials = block.trial;
                    
                    disp('Loading lick data...');
                    try
                        L=load(dat.expFilePath(obj.expRef, 'Timeline', 'm'));
                        tseries = L.Timeline.rawDAQTimestamps;
                        lickseries = L.Timeline.rawDAQData(:,7);
                        
                        for t=1:block.numCompletedTrials
                            start = trials(t).trialStartedTime;
                            finish = trials(t).trialEndedTime;
                            idx = (start < tseries) & (tseries < finish);
                            
                            obj.data.lickenergy(t,1) = sum(lickseries(idx).^2);
                        end
                    catch
                        warning('No lick data found.. Setting to NaNs');
                        obj.data.lickenergy = nan(block.numCompletedTrials,1);
                    end
                    
                case 'pupil'
                    block = dat.loadBlock(obj.expRef);
                    trials = block.trial;
                    
                    disp('Loading eye movie...');
                    
                    try
                        vidFilename = ['\\zserver\Data\EyeCamera\' obj.expRef(14:end) '\' obj.expRef(1:10) '\' obj.expRef(12) '\eye.mj2'];
                        v = VideoReader(vidFilename);
                        
                        for t=1:block.numCompletedTrials
                            %                             disp(t);
                            start = trials(t).trialStartedTime;
                            finish = trials(t).trialEndedTime;
                            
                            v.CurrentTime = start; a=0; pixCounts=0;
                            while v.CurrentTime < finish
                                frame = 255 - (v.readFrame*5.5);
                                frame(frame>0)=1;
                                frame = frame(50:110,80:180);
                                pixCounts = pixCounts + sqrt(sum(frame(:)));
                                a=a+1;
                            end
                            
                            obj.data.pupilenergy(t,1) = pixCounts/a;
                        end
                    catch
                        warning('No eye data found.. Setting to NaNs');
                        obj.data.pupilenergy = nan(block.numCompletedTrials,1);
                        
                    end
            end
            
            disp('Done!');
            
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
        
        function paramCovSimulate(obj)
            %Plots covariance/variance of prms as a function of
            %varying stimulus and behavioural configurations. This is to
            %explore other ways of setting up the task such that the
            %prm estimation can be less variable/less covariable.
            if isempty(obj.prmFits)
                error('Fit first');
            end
            
            %1) Same prms but vary the distribution of contrasts
            if obj.ContrastDimensions == 1
                n = obj.prmFits(end-1);
                c50 = obj.prmFits(end);
                cfn = @(c)(c.^n)./(c.^n + c50.^n);
                numTrials = length(obj.data.response);
                
                widths = linspace(0.01,0.5,50);
                
                meanVarS = []; cv = [];
                for w = 1:length(widths)
                    varS = [];
                    for iter = 1:50
                        crnd = widths(w)*randn(numTrials,1);
                        c = zeros(length(crnd),2);
                        c(crnd<0,1) = -crnd(crnd<0);
                        c(crnd>0,2) = crnd(crnd>0);
                        
                        H = obj.hessian(obj.prmFits,c);
                        
                        F = -sum(H,3);
                        cov = inv(F);
                        varS(iter,:) = [cov(1,1) cov(2,2) cov(3,3) cov(4,4) cov(3,1) cov(2,1) cov(4,3)];
                        
                    end
                    meanVarS(w,:) = mean(varS);
                    
                    corrSLBL(w,1) = meanVarS(w,5)/sqrt(meanVarS(w,1)*meanVarS(w,3));
                    
                    cv(w,1) = sqrt(meanVarS(w,1))/obj.prmFits(1);
                    cv(w,2) = sqrt(meanVarS(w,2))/obj.prmFits(3);
                    cv(w,3) = sqrt(meanVarS(w,3))/obj.prmFits(2);
                    cv(w,4) = sqrt(meanVarS(w,4))/obj.prmFits(4);
                end
                figure;
                subplot(3,1,1);
                plot(widths,[meanVarS corrSLBL],'LineWidth',2); xlabel('Standard Deviation of Contrasts');
                legend({'Var[B_L]','Var[B_R]','Var[S_L]','Var[S_R]','Cov[B_L,S_L]','Cov[B_L,B_R]','Cov[S_L,S_R]','Corr[B_L,S_L]'});
                
                hold on;
                this = std(diff(obj.data.stimulus,[],2));
                line([this,this],[-1 1]);
            end
            
            %2) Same contrasts but vary GovsNoGo bias
            fitted_p = obj.prmFits;
            Biases = linspace(-5,5,50);
            varS = []; CVs = [];
            for bias = 1:length(Biases)
                new_p = fitted_p;
                new_p([1 3]) = Biases(bias);
                H = obj.hessian(new_p,obj.data.stimulus);
                F = -sum(H,3);
                cov = inv(F);
                varS(bias,:) = [cov(1,1) cov(2,2) cov(3,3) cov(4,4) cov(3,1) cov(2,1) cov(4,3)];
                CVs(bias,1) = sqrt(cov(1,1))/new_p(1);
                CVs(bias,2) = sqrt(cov(2,2))/new_p(3);
                CVs(bias,3) = sqrt(cov(3,3))/new_p(2);
                CVs(bias,4) = sqrt(cov(4,4))/new_p(4);
            end
            corrS = varS(:,5)./sqrt(varS(:,1).*varS(:,3));
            subplot(3,1,2);
            pGO = 1- 1./(1+2*exp(Biases)); xlim([0 1]);
            plot(pGO,[varS,corrS],'LineWidth',2); xlabel('pGO');
            %             legend({'cv[B_L]','cv[B_R]','cv[S_L]','cv[S_R]'});
            
            this = mean(fitted_p([1 3]));
            this = 1- 1./(1+2*exp(this));
            line([this,this],[-1 1]);
            
            %3) Same contrast but vary LvNG bias
            fitted_p = obj.prmFits;
            Biases = linspace(-5,5,50);
            varS = [];
            for bias = 1:length(Biases)
                new_p = fitted_p;
                new_p(1) = Biases(bias);
                H = obj.hessian(new_p,obj.data.stimulus);
                F = -sum(H,3);
                cov = inv(F);
                varS(bias,:) = [cov(1,1) cov(2,2) cov(3,3) cov(4,4) cov(3,1) cov(2,1) cov(4,3)];
            end
            corrS = varS(:,5)./sqrt(varS(:,1).*varS(:,3));
            subplot(3,1,3);
            plot(Biases,[varS corrS],'LineWidth',2); xlabel(['bL . bR=' num2str(fitted_p(3))]);
            this = fitted_p(1);
            %             this = 1- 1./(1+2*exp(this));
            line([this,this],[-1 1]);
        end
        
        function plotPspace_mnrvsnested(obj)
            %use mnrfit to fit rudimentary models MNR and NESTED and
            %compare directly
            cont = obj.data.stimulus;
            resp = obj.data.response;
            resp_nes = resp; resp_nes(resp_nes==3)=0; resp_nes=resp_nes+1;
            
            B_mnr = mnrfit(cont,resp);
            B_nes = mnrfit(cont,resp_nes,'model','hierarchical');
            
            cVals = linspace(0,0.54,1000);
            [cr,cl]=meshgrid(cVals);
            cont = [cl(:) cr(:)];
            
            P_mnr = mnrval(B_mnr,cont);
            P_nes = mnrval(B_nes,cont,'model','hierarchical');
            P_nes = [P_nes(:,2:3) P_nes(:,1)];
            
            %Plot
            figure('color','w');
            labels = {'pL','pR','pNG'};
            %             cVals=log(cVals);
            for r = 1:3
                subplot(1,3,r);
                [~,ax]=contour(cVals,cVals,reshape(P_mnr(:,r),size(cl)));
                ax.LineWidth=1; hold on;
                [~,ax]=contour(cVals,cVals,reshape(P_nes(:,r),size(cl)));
                ax.LineWidth=2;
                ax.LineStyle=':';
                
                set(gca,'ydir','normal','box','off');
                xlabel('CR'); ylabel('CL'); title(labels{r}); axis square; caxis([0 1]);
            end
            %             keyboard;
        end
        
        function bootstrapFit(obj)
            if isempty(obj.ZL); error('Set model first'); end
            
            numIter = 300;
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
            bootstrap_params = nan(numIter,length(obj.prmLabels));
            
            T = struct2table(obj.data);
            for iter = 1:numIter
                datTemp = datasample(T,height(T));
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(datTemp), datTemp.response));
                bootstrap_params(iter,:) = fmincon(objective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
            end
            
            gplotmatrix(bootstrap_params,[],[],[],[],[],[],[],obj.prmLabels);
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
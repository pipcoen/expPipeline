function [logOddsLR, logOddsTO] = GLMMultiModels(obj, tag, P)
if ~exist('P', 'var'); obj.modelString = tag; end
% obj.blockData = prc.combineBlocks(obj.blockData, obj.blockData.selectedTrials);
[obj.blockData.audValues, uniA] = deal(unique(obj.blockData.audDiff));
[obj.blockData.visValues, uniV] = deal(unique(obj.blockData.visDiff));
uniA = sort(uniA, 'descend'); uniV = sort(uniV, 'descend');

[audGrid,visGrid] = meshgrid(uniA,uniV);
comb = unique([obj.blockData.visDiff obj.blockData.audDiff], 'rows');
switch tag
    case 'eval'; visDiff = obj.evalPoints(:,1); audDiff = obj.evalPoints(:,2);
    otherwise, audDiff = obj.blockData.audDiff; visDiff = obj.blockData.visDiff;
end
if contains(obj.modelString, 'Nest'); nested = 1; else; nested = 0; end
if contains(obj.modelString, 'Conf'); addConf = 1; else; addConf = 0; end
confTrials = (sign(audDiff).*sign(visDiff))==-1;
repAud = repmat({audDiff},1,length(uniA));
repVis = repmat({visDiff},1,length(uniV));
audTags = arrayfun(@(x) [num2str(x) 'Aud'], uniA, 'uni', 0);

getC50 = @(visData,N,C50) (visData.^N)./(visData.^N + C50^N);

switch obj.modelString        
    case {'biasOnly';'biasOnlyNest';'VisOnly';'AudOnly';'AudDom'; 'SimpLog';'SimpLogNest'; 'SimpLogNestConf'; 'SimpLogBiasTONest'}
        if contains(obj.modelString, 'VisOnly'); vOnly = 0; else; vOnly = 1; end
        if contains(obj.modelString, 'AudOnly'); aOnly = 0; else; aOnly = 1; end
        if contains(obj.modelString, 'biasOnly'); bOnly = 0; else; bOnly = 1; end
        if any(contains(obj.modelString, {'BiasTO'; 'BiasTO'})); bOnlyTO = 0; else; bOnlyTO = 1; end
        if contains(obj.modelString, 'AudDom'); domIdx = confTrials; else; domIdx = zeros(length(confTrials),1); end
        obj.prmLabels = [{'bias';'visScale';'N';'C50'}; audTags];
        
        TOSrt = length(obj.prmLabels);
        if nested; obj.prmLabels = [obj.prmLabels; cellfun(@(x) [x 'TO'], obj.prmLabels(~contains(obj.prmLabels, {'N';'C50'})), 'uni', 0)]; end
        if addConf; obj.prmLabels = [obj.prmLabels; 'confWeight']; end
        if exist('P', 'var')
            visContributionLR = (P(2)*(getC50(abs(visDiff), P(3), P(4))).*sign(visDiff)).*(~domIdx);
            audContributionLR = arrayfun(@(x,y,z) x*(y{1}==z), P(5:TOSrt), repAud, uniA', 'uni', 0);
            logOddsLR = P(1)+visContributionLR*aOnly*bOnly + sum(cell2mat(audContributionLR),2)*vOnly*bOnly;
            
            if nested
                visContributionTO = P(TOSrt+2)*(getC50(abs(visDiff), P(3), P(4)));
                audContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P(TOSrt+3:TOSrt+2+length(uniA)), repAud, uniA', 'uni', 0);
                logOddsTO = P(TOSrt+1)+visContributionTO*bOnly*bOnlyTO +  sum(cell2mat(audContributionTO),2)*bOnly*bOnlyTO;
                if addConf; logOddsTO = logOddsTO+P(end)*double(confTrials); end
            end
            
        end
        obj.evalPoints = [repmat(linspace(-max(abs(uniV)),max(abs(uniV)),200)', length(uniA),1), reshape(repmat(uniA,1,200)',600,1)];
        obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
        
    case {'SimpLogSplit'}
        obj.prmLabels = [{'bias';'visScaleR';'visScaleL';'N';'C50'}; audTags];
        TOSrt = length(obj.prmLabels);
        if nested; obj.prmLabels = [obj.prmLabels, cellfun(@(x) [x 'TO'], obj.prmLabels, 'uni', 0), 'conflictTO']; end
        if exist('P', 'var')
            freeP = ones(8,1);
            P = P+obj.prmInit;
            visContributionLR = (freeP(2)*P(2)*(getC50(abs(visDiff.*(visDiff>0)), freeP(4)*P(4), freeP(5)*P(5))) + ...
                freeP(3)*P(3)*getC50(abs(visDiff.*(visDiff<0)), freeP(4)*P(4), freeP(5)*P(5)));
            audContributionLR = arrayfun(@(x,y,z) x* (y{1}==z), P(6:end).*freeP(6:end)', repAud, uniA', 'uni', 0);
            logOddsLR = freeP(1)*P(1)+visContributionLR+sum(cell2mat(audContributionLR),2);
            if nested
                visContributionTO = P(TOSrt+2)*(sqrt(abs(visDiff).*(visDiff>0))) + P(TOSrt+3)*(sqrt(abs(visDiff).*(visDiff<0)));
                audContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P(TOSrt+4:TOSrt+TOSrt), repAud, uniA', 'uni', 0);
                logOddsTO = P(TOSrt+1)+visContributionTO+sum(cell2mat(audContributionTO),2) + P(end)*double(confTrials);
            end
        end
        obj.evalPoints = [repmat(linspace(-max(abs(uniV)),max(abs(uniV)),200)', length(uniA),1), reshape(repmat(uniA,1,200)',600,1)];
        
    case {'SimpLogSplitDelta1';'SimpLogSplitDelta2';'SimpLogSplitDelta3';'SimpLogSplitDelta4';'SimpLogSplitDelta5';...
            'SimpLogSplitDelta6';'SimpLogSplitDelta7';'SimpLogSplitDelta8'}
        obj.prmLabels = {'delta'};
        pInit = obj.prmInit;
        freeP = str2double(obj.modelString(end));
        TOSrt = length(obj.prmLabels);
        if exist('P', 'var')
            visContributionLR = ((pInit(2)+(P(1)*freeP==2))*(getC50(abs(visDiff.*(visDiff>0)), (pInit(4)+(P(1)*freeP==4)), (pInit(5)+(P(1)*freeP==5)))) + ...
               (pInit(3)+(P(1)*freeP==3))*getC50(abs(visDiff.*(visDiff<0)), (pInit(4)+(P(1)*freeP==4)), (pInit(5)+(P(1)*freeP==5))));
            audContributionLR = arrayfun(@(x,y,z) x*(y{1}==z), [(pInit(6)+(P(1)*freeP==6)), (pInit(7)+(P(1)*freeP==7)), (pInit(8)+(P(1)*freeP==8))], repAud, uniA', 'uni', 0);
            logOddsLR = (pInit(1)+(P(1)*freeP==1))+visContributionLR+sum(cell2mat(audContributionLR),2);
            if nested
                visContributionTO = P(TOSrt+2)*(sqrt(abs(visDiff).*(visDiff>0))) + P(TOSrt+3)*(sqrt(abs(visDiff).*(visDiff<0)));
                audContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P(TOSrt+4:TOSrt+TOSrt), repAud, uniA', 'uni', 0);
                logOddsTO = P(TOSrt+1)+visContributionTO+sum(cell2mat(audContributionTO),2) + P(end)*double(confTrials);
            end
        end
        obj.evalPoints = [repmat(linspace(-max(abs(uniV)),max(abs(uniV)),200)', length(uniA),1), reshape(repmat(uniA,1,200)',600,1)];
        
    case {'SimpEmp'; 'SimpEmpNest'; 'SimpEmpNestConf'}
        obj.prmLabels = ['bias'; arrayfun(@(x) [num2str(x) 'Vis'], uniV, 'uni', 0); audTags];
        TOSrt = length(obj.prmLabels);
        if nested; obj.prmLabels = [obj.prmLabels; cellfun(@(x) [x 'TO'], obj.prmLabels, 'uni', 0)]; end
        if contains(obj.modelString, 'Conf'); obj.prmLabels = [obj.prmLabels; 'ConfIdent']; end
        if exist('P', 'var')
            visContributionLR = arrayfun(@(x,y,z) x*(y{1}==z), P(2:(length(uniV)+1)), repVis, uniV', 'uni', 0);
            audContributionLR = arrayfun(@(x,y,z) x*(y{1}==z), P((length(uniV)+2):TOSrt), repAud, uniA', 'uni', 0);
            logOddsLR = P(1)+sum(cell2mat(visContributionLR), 2)+sum(cell2mat(audContributionLR),2);
            
            if nested
                visContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P(2+TOSrt:(length(uniV)+1)+TOSrt), repVis, uniV', 'uni', 0);
                audContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P((length(uniV)+2)+TOSrt:TOSrt+TOSrt), repAud, uniA', 'uni', 0);
                logOddsTO = P(TOSrt+1)+sum(cell2mat(visContributionTO), 2)+sum(cell2mat(audContributionTO),2);
                if contains(obj.modelString, 'Conf'); logOddsTO = logOddsTO + P(end)*double(confTrials); end
            end
        end
        obj.evalPoints = [visGrid(:) audGrid(:)];

        
    case {'FullEmp'; 'FullEmpNest'}
        obj.prmLabels = ['bias'; cellfun(@(x) [sprintf('%0.1f', x), 'VisAud'], num2cell(comb,2), 'uni', 0)];
        TOSrt = length(obj.prmLabels);
        if nested; obj.prmLabels = [obj.prmLabels; cellfun(@(x) [x 'TO'], obj.prmLabels, 'uni', 0)]; end
        if exist('P', 'var')
            repeatedValues = repmat({[visDiff, audDiff]},1,size(comb,1));
            stimulusContributionsLR = arrayfun(@(x,y,z) x.*(all(y{1}==z{1},2)),P(2:TOSrt),repeatedValues, num2cell(comb,2)', 'uni', 0);
            logOddsLR = P(1)+sum(cell2mat(stimulusContributionsLR), 2);
            
            if nested
                stimulusContributionsTO = arrayfun(@(x,y,z) x.*(all(y{1}==z{1},2)),P(2+TOSrt:TOSrt+TOSrt),repeatedValues, num2cell(comb,2)', 'uni', 0);
                logOddsTO = P(1+TOSrt)+sum(cell2mat(stimulusContributionsTO), 2);
            end
        end
        obj.evalPoints = comb;
        
    otherwise, error('modelString not recgonized');
end
if isempty(obj.prmBounds) || size(obj.prmBounds,2)~= length(obj.prmLabels)
    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
end
if any(strcmp(obj.prmLabels, 'N')); obj.prmBounds(:, strcmp(obj.prmLabels, 'N')) = [0;3]; end
if any(strcmp(obj.prmLabels, 'C50')); obj.prmBounds(:, strcmp(obj.prmLabels, 'C50')) = [0.001;0.9]; end



if isempty(obj.prmInit)
    obj.prmInit = zeros(1,size(obj.prmBounds,2));
end

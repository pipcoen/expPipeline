function [logOddsLR, logOddsTO] = GLMMultiModels(obj, tag, P)
if ~exist('P', 'var'); obj.modelString = tag; end
% obj.blockData = prc.combineBlocks(obj.blockData, obj.blockData.selectedTrials);
[obj.blockData.audValues, uniA] = deal(unique(obj.blockData.audDiff));
[obj.blockData.visValues, uniV] = deal(unique(obj.blockData.visDiff));

[audGrid,visGrid] = meshgrid(uniA,uniV);
comb = unique([obj.blockData.visDiff obj.blockData.audDiff], 'rows');
switch tag
    case 'eval'; visDiff = obj.evalPoints(:,1); audDiff = obj.evalPoints(:,2);
    otherwise audDiff = obj.blockData.audDiff; visDiff = obj.blockData.visDiff;
end
if contains(obj.modelString, 'Nest'); nested = 1; else; nested = 0; end
confTrials = (sign(audDiff).*sign(visDiff))==-1;

switch obj.modelString
    case 'Simp-log'
        obj.prmLabels = [{'bias','visScale'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(sort(uniA, 'descend')), 'uni', 0)'];
        if exist('P', 'var')
            visContribution = P(2)*(sqrt(abs(visDiff)).*sign(visDiff));
            audContribution = arrayfun(@(x,y,z) x*(y{1}==z), P(3:end), repmat({audDiff},1,length(uniA)), sort(uniA, 'descend')', 'uni', 0);
            logOddsLR = P(1)+visContribution+sum(cell2mat(audContribution),2);
        end
        obj.evalPoints = [repmat(linspace(-max(abs(uniV)),max(abs(uniV)),200)', length(uniA),1), reshape(repmat(uniA,1,200)',600,1)];
        
    case {'SqrtLogisticSplit'; 'SqrtLogisticSplitNest'}
        obj.prmLabels = [{'bias','visScaleR','visScaleL'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(sort(uniA, 'descend')), 'uni', 0)'];
        TOSrt = length(obj.prmLabels);
        if nested; obj.prmLabels = [obj.prmLabels, cellfun(@(x) [x 'TO'], obj.prmLabels, 'uni', 0), 'conflictTO']; end
        if exist('P', 'var')
            visContributionLR = P(2)*(sqrt(abs(visDiff).*(visDiff>0))) + P(3)*(sqrt(abs(visDiff).*(visDiff<0)));
            audContributionLR = arrayfun(@(x,y,z) x*(y{1}==z), P(4:TOSrt), repmat({audDiff},1,length(uniA)), sort(uniA, 'descend')', 'uni', 0);
            logOddsLR = P(1)+visContributionLR+sum(cell2mat(audContributionLR),2);
            
            if nested
                visContributionTO = P(TOSrt+2)*(sqrt(abs(visDiff).*(visDiff>0))) + P(TOSrt+3)*(sqrt(abs(visDiff).*(visDiff<0)));
                audContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P(TOSrt+4:TOSrt+TOSrt), repmat({audDiff},1,length(uniA)), sort(uniA, 'descend')', 'uni', 0);
                logOddsTO = P(TOSrt+1)+visContributionTO+sum(cell2mat(audContributionTO),2) + P(end)*double(confTrials);
            end
        end
        obj.evalPoints = [repmat(linspace(-max(abs(uniV)),max(abs(uniV)),200)', length(uniA),1), reshape(repmat(uniA,1,200)',600,1)];
        
    case {'SimpEmp'; 'SimpEmpNest'}
        obj.prmLabels = ['bias'; cellfun(@(x) [num2str(x) 'Vis'], num2cell(uniV), 'uni', 0); cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)];
        TOSrt = length(obj.prmLabels);
        if nested; obj.prmLabels = [obj.prmLabels; cellfun(@(x) [x 'TO'], obj.prmLabels, 'uni', 0)]; end
        if exist('P', 'var')
            visContributionLR = arrayfun(@(x,y,z) x*(y{1}==z), P(2:(length(uniV)+1)), repmat({visDiff},1,length(uniV)), sort(uniV, 'descend')', 'uni', 0);
            audContributionLR = arrayfun(@(x,y,z) x*(y{1}==z), P((length(uniV)+2):TOSrt), repmat({audDiff},1,length(uniA)), sort(uniA, 'descend')', 'uni', 0);
            logOddsLR = P(1)+sum(cell2mat(visContributionLR), 2)+sum(cell2mat(audContributionLR),2);
            
            if nested
                visContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P(2+TOSrt:(length(uniV)+1)+TOSrt), repmat({visDiff},1,length(uniV)), sort(uniV, 'descend')', 'uni', 0);
                audContributionTO = arrayfun(@(x,y,z) x*(y{1}==z), P((length(uniV)+2)+TOSrt:TOSrt+TOSrt), repmat({audDiff},1,length(uniA)), sort(uniA, 'descend')', 'uni', 0);
                logOddsTO = P(TOSrt+1)+sum(cell2mat(visContributionTO), 2)+sum(cell2mat(audContributionTO),2);
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
end
if isempty(obj.prmBounds) || size(obj.prmBounds,2)~= length(obj.prmLabels)
    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
end
if isempty(obj.prmInit)
    obj.prmInit = zeros(1,size(obj.prmBounds,2));
end

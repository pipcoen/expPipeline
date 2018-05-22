function logLik = GLMMultiModels(obj, tag, P)
if ~exist('P', 'var'); obj.modelString = tag; end
[obj.blockData.audValues, uniA] = deal(unique(obj.blockData.audDiff));
[obj.blockData.visValues, uniV] = deal(unique(obj.blockData.visDiff));
obj.evalPoints = [repmat(linspace(-max(abs(uniV)),max(abs(uniV)),200)', length(uniA),1), reshape(repmat(uniA,1,200)',600,1)];

% [audGrid,visGrid] = meshgrid(uniA,uniV);
if strcmpi(tag, 'eval'); visDiff = obj.evalPoints(:,1); audDiff = obj.evalPoints(:,2);
else, audDiff = obj.blockData.audDiff; visDiff = obj.blockData.visDiff;
end

switch obj.modelString
    case 'SqrtLogisticSplitDelta'
        obj.prmLabels = [{'bias','visScaleR','visScaleL'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(sort(uniA, 'descend')), 'uni', 0)'];
        obj.prmLabels = cellfun(@(x) [x 'Delta'], obj.prmLabels, 'uni', 0);
        if exist('P', 'var')
            P = P+obj.prmInit;
            visContribution = P(2)*(sqrt(abs(visDiff).*(visDiff>0))) + P(3)*(sqrt(abs(visDiff).*(visDiff<0)));
            audContribution = arrayfun(@(x,y,z) x*(y{1}==z), P(4:end), repmat({audDiff},1,length(uniA)), sort(uniA, 'descend')', 'uni', 0);
            logLik = P(1)+visContribution+sum(cell2mat(audContribution),2);
        end
end
if isempty(obj.prmBounds) || size(obj.prmBounds,2)~= length(obj.prmLabels)
    obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
end
if isempty(obj.prmInit)
    obj.prmInit = zeros(1,size(obj.prmBounds,2));
end



% switch
% 
% uniV = blk.visValues;
% comb = unique([obj.blockData.visDiff obj.blockData.audDiff], 'rows');
% numA = length(uniA);
% numV = length(uniV);
% numC = length(comb);
% 
% 
% maxContrast = max(abs(uniV));
% audRepMat = repmat(uniA,1,200)';
% obj.evalPoints = [repmat(linspace(-maxContrast,maxContrast,200)', numA,1), audRepMat(:)];
% [aGrid,vGrid] = meshgrid(uniA,uniV);
% 'SqrtLogisticSplitDelta'
% obj.prmLabels = [{'bias','visScaleR','visScaleL'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)'];
% obj.prmLabels = cellfun(@(x) [x 'Delta'], obj.prmLabels, 'uni', 0);
% obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
% obj.inputFun = @(b)({b.visDiff; b.audDiff; b.prmInit});
% obj.modelFun = @(P,in) (P(1)+in{3}(1) + (P(2)+in{3}(2))*sqrtDiff(in{1}.*(in{1}>0)) + (P(3)+in{3}(3))*sqrtDiff(in{1}.*(in{1}<0)) + ...
%     indiSum(P(4:(3+numA))+in{3}(4:(3+numA)),in{2}));
% 
% 

%                 case 'SimpleLogistic'
%                     obj.prmLabels = {'bias','visScale','audScale'};
%                     obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
%                     obj.inputFun = @(b)([b.visDiff, b.audDiff]);
%                     obj.modelFun = @(P,in) (P(1)+ P(2)*in(:,1) + P(3)*in(:,2));
%                 case 'SqrtLogistic'
%                     obj.prmLabels = [{'bias','visScale'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)'];
%                     obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
%                     obj.inputFun = @(b)([b.visDiff, b.audDiff]);
%                     obj.modelFun = @(P,in) (P(1)+ P(2)*sqrtDiff(in(:,1))+indiSum(P(3:(2+numA)),in(:,2)));
%                 case 'SqrtLogisticSplit'
%                     obj.prmLabels = [{'bias','visScaleR','visScaleL'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)'];
%                     obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
%                     obj.inputFun = @(b)([b.visDiff, b.audDiff]);
%                     obj.modelFun = @(P,in) (P(1)+ P(2)*sqrtDiff(in(:,1).*(in(:,1)>0))+P(3)*sqrtDiff(in(:,1).*(in(:,1)<0))+indiSum(P(4:(3+numA)),in(:,2)));
%                 case 'SqrtLogisticSplit'
%                     obj.prmLabels = [{'bias','visScaleR','visScaleL'},cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)'];
%                     obj.prmLabels = cellfun(@(x) [x 'Delta'], obj.prmLabels, 'uni', 0);
                    
%                     obj.modelFun = @(P,in) (P(1)+in{3}(1) + (P(2)+in{3}(2))*sqrtDiff(in{1}.*(in{1}>0)) + (P(3)+in{3}(3))*sqrtDiff(in{1}.*(in{1}<0)) + ...
%                         indiSum(P(4:(3+numA))+in{3}(4:(3+numA)),in{2}));
%                 case 'Simp-emp'
%                     allValues = [cellfun(@(x) [num2str(x) 'Vis'], num2cell(uniV), 'uni', 0); cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0)];
%                     obj.prmLabels = ['bias'; allValues(:)]';
%                     obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
%                     obj.inputFun = @(b)([b.visDiff, b.audDiff]);
%                     obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numV)),in(:,1)) + indiSum(P((2+numV):(1+numV+numA)),in(:,2)));
%                     obj.evalPoints = [vGrid(:) aGrid(:)];
%                 case 'Simp-aud'
%                     allValues = cellfun(@(x) [num2str(x) 'Aud'], num2cell(uniA), 'uni', 0);
%                     obj.prmLabels = ['bias'; allValues(:)]';
%                     obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
%                     obj.inputFun = @(b)([b.visDiff, b.audDiff]);
%                     obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numA)),in(:,2)));
%                     obj.evalPoints = [vGrid(:) aGrid(:)];
%                 case 'Simp-vis'
%                     allValues = cellfun(@(x) [num2str(x) 'Vis'], num2cell(uniV), 'uni', 0);
%                     obj.prmLabels = ['bias'; allValues(:)]';
%                     obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
%                     obj.inputFun = @(b)([b.visDiff, b.audDiff]);
%                     obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numV)),in(:,1)));
%                     obj.evalPoints = [vGrid(:) aGrid(:)];
%                 case 'Full-emp'
%                     allValues = [cellfun(@(x) [sprintf('%0.1f', x), 'VisAud'], num2cell(comb,2), 'uni', 0)];
%                     obj.prmLabels = ['bias'; allValues]';
%                     obj.prmBounds = repmat([-inf; inf], 1, length(obj.prmLabels));
%                     obj.inputFun = @(b)([b.visDiff, b.audDiff]);
%                     obj.modelFun = @(P,in)(P(1) + indiSum(P(2:(1+numC)),in));
%                     obj.evalPoints = comb;


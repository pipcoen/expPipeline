classdef passiveData
    
    properties
        mNam;
        rDat;
        stEn;
        sSrt;
        spkT;
        spkN;
        spkA;
        cPln;
        avIx;
    end
    
    methods
        function obj = passiveData

        end
        
        function p = tempRCells(obj)
            for i = 1:length(obj.mNam)
                for j = 1:length(obj.spkN{i})
                    cohS = double([squeeze(obj.spkN{i}(j,1,1,:)); ...
                        squeeze(obj.spkN{i}(j,2,2,:))]);
                    incS = double([squeeze(obj.spkN{i}(j,2,1,:)); ...
                        squeeze(obj.spkN{i}(j,1,2,:))]);
                    p{i}(j,1) = anova1([cohS incS], [], 'off');
                end
            end
        end
        
        function [p, spkN] = visRCells(obj)
            for i = 1:length(obj.mNam)
                for j = 1:length(obj.spkN{i})
                    vStm = double(squeeze(obj.spkN{i}(j,:,1:4,:)));
                    nStm = double(squeeze(obj.spkN{i}(j,:,5,:)));
                    [~, p{i}(j,1)] = ttest2(vStm(:), nStm(:));
                end
                spkN{i} = obj.spkN{i}(p{i}*length(p{i})<0.01,:,:,:);
            end
        end
    end
end
function collectOriTuning
fLst = rdir([savePath('') '*\**\5\*flu.mat']);
oriT = struct;
for i = 1:length(fLst)
    load(fLst(i).name);
    load(savePath('problock', flu.mNam, flu.rDat, [], '5'));
    oriT(i).mNam = flu.mNam;
    oriT(i).rDat = flu.rDat;
    oriT(i).stEn = blk.stEn;
    oriT(i).sPre = blk.sPre;
    oriT(i).oriA = blk.oriA;
    oriT(i).cPln = flu.cPln;
    
    spkT = cell(length(flu.cPln),1);
    spkA = cell(length(flu.cPln),1);
    for j = 1:length(flu.cPln)
        tVal = flu.frmT(flu.cPln(j),flu.spkT{j})';
        [spkT{j}, sIdx] = sort(tVal); %#ok<TRSRT>
        spkA{j} = flu.spkA{j,1}(sIdx)';
    end
    obj.spkT{i,1} = spkT;
    obj.spkA{i,1} = spkA;
end
end

%         function [cFit] = oriStrength(obj)
%             for i = 1:length(obj.mNam)
%                 sPre = obj.sPre{i};
%                 oriA = obj.oriA{i};
%                 oris = sort(unique(obj.oriA{i}));
%                 for j = 1:length(oris)
%                     idxT = sPre(oriA==oris(j),:);
%                     spkR(:,j) = mean(normStimSpks(idxT, obj.spkT{i}, 0.5),2);
%                 end
%                 for j = 1:size(spkR,1)
%                     spkR(j,:) = spkR(j,:)-min(spkR(j,:));
%                     dirS(j,1) = 1-abs(sum(spkR(j,:)'.*exp((2*1i*pi/360)*oris)))./sum(spkR(j,:));
%                     oriS(j,1) = 1-abs(sum(spkR(j,:)'.*exp((4*1i*pi/360)*oris)))./sum(spkR(j,:));
%                 end
%             end
%         end
%
%     end
% end
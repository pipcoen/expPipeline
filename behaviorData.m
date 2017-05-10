classdef behaviorData
    
    properties
        mNam;
        rDat;
        fCor;
        tNum;
    end
    
    methods
        function obj = behaviorData
           fLst = rdir([savePath('') '*\**\*blk.mat']);
           idx = 0;
           obj.rDat = '';
           for i = 1:length(fLst)
               load(fLst(i).name);
               if ~any([1 2 4] == str2double(blk.eIdx)); continue; end
               if length(blk.fBck) < 50
                   fprintf('Skipping %s %s due to low trial number\n', ...
                       blk.mNam, blk.rDat);
                   continue;
               end
               
               idx = idx+1;
               obj.mNam{idx,1} = blk.mNam;
               obj.rDat(idx,:) = blk.rDat;
               obj.fCor(idx,1) = (mean(blk.fBck)+1)/2;
               obj.tNum(idx,1) = length(blk.fBck);
           end
        end
        
        function fracCor(obj, mNam)
            perf = obj.fCor(strcmp(obj.mNam, mNam));
            date = obj.rDat(strcmp(obj.mNam, mNam),:);
            nDay = daysact(date(1,:), date);
            plot(nDay, perf, '.-');
        end
    end
end
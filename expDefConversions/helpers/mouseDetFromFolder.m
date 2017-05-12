function [subject, expDate, sessionNum] = mouseDetFromFolder(foldrNam)
sepPoint = find(foldrNam =='/' | foldrNam =='\');
sepPoint(sepPoint < strfind(foldrNam, 'expInfo')) = [];
subject = foldrNam(sepPoint(1)+1:sepPoint(2)-1);
expDate = foldrNam(sepPoint(2)+1:sepPoint(3)-1);
sessionNum = foldrNam(sepPoint(3)+1:end);
end
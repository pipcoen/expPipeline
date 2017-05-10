function [mouseNam, exptDate, sessnNum] = mouseDetFromFolder(foldrNam)
sepPoint = find(foldrNam =='/' | foldrNam =='\');
sepPoint(sepPoint < strfind(foldrNam, 'expInfo')) = [];
mouseNam = foldrNam(sepPoint(1)+1:sepPoint(2)-1);
exptDate = foldrNam(sepPoint(2)+1:sepPoint(3)-1);
sessnNum = foldrNam(sepPoint(3)+1:end);
end
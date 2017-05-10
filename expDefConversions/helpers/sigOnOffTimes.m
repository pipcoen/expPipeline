function [onOff] = sigOnOffTimes(tVal, tTim)
if tVal(end)==1; tVal(end) = []; end
onOff = [tTim(tVal==1)' tTim(strfind(tVal, [1 0])+1)'];
end
function stdPatch(xDat, mnVal, stdVal, patchColor, alphaValue, edgeOnly)
if ~exist('alphaValue', 'var');alphaValue = 0.5; end
if ~exist('edgeOnly', 'var');edgeOnly = 0; end
nVal = numel(xDat);
mnVal = mnVal(:)';
stdVal = stdVal(:)';
patchIdx = ([1:nVal nVal:-1:1 1]);
if ~edgeOnly
    patch(xDat(patchIdx), [mnVal+stdVal mnVal(end:-1:1)-stdVal(end:-1:1) mnVal(1)+stdVal(1)], patchColor, 'FaceAlpha',alphaValue, 'EdgeColor', 'none')
else, plot(xDat(patchIdx), [mnVal+stdVal mnVal(end:-1:1)-stdVal(end:-1:1) mnVal(1)+stdVal(1)], 'color', patchColor)
end
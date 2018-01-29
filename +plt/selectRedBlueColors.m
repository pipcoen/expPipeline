function colorChoices = selectRedBlueColors(inputValues, zeroTag)
if ~exist('zeroTag', 'var'); zeroTag = 0.5*ones(1,3); end
allColors = plt.redblue(255);
maxLength = max(abs(inputValues));
colorChoices = zeros(length(inputValues),3);
fractionalPosition = inputValues./maxLength;
colorChoices(inputValues>0,:) = allColors(128-(round(fractionalPosition(inputValues>0)*127)),:);
colorChoices(inputValues<0,:) = allColors(128+(round(fractionalPosition(inputValues<0)*-127)),:);
if any(inputValues == 0); colorChoices(inputValues==0,:) = zeroTag; end
colorChoices = flip(colorChoices, 1);
end
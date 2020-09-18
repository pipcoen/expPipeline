function blk = getDefinedSubset(blk, subsetTag, minMaxVis)
%% Function to filter the trials in a block based on some predefined tags

%INPUTS(default values)
%blk(required)---------The struct array of block files to be filtered
%subsetTag(required)---A tag to indicate which subset of trials to select
%   'blnk'-------------------blank trials (where aud is in the center and visual contrast is zero)
%   'vl'---------------------visual trials (aud in the center) where vis initial azimuth in on the left
%   'vr'---------------------visual trials (aud in the center) where vis initial azimuth in on the right
%   'al'---------------------auditory trials (visual contrast is zero) where aud initial azimuth in on the left
%   'ar'---------------------auditory trials (visual contrast is zero) where aud initial azimuth in on the right
%   'cohl'-------------------coherent trials (both stimuli informative) where initial azimuth in on the left
%   'cohr'-------------------coherent trials (both stimuli informative) where initial azimuth in on the right
%   'conl'-------------------conflict trials (both stimuli informative) where aud initial vis/aud azimuth in on the left/right
%   'conr'-------------------conflict trials (both stimuli informative) where aud initial aud/vis azimuth in on the right/left

%OUTPUTS
%blk-------------Filtered block

if ~exist('minMaxVis', 'var'); minMaxVis = [0 1]; end
if ~any(strcmpi(subsetTag, {'blnk';'al';'ar'}))
    blk = prc.filtBlock(blk, blk.tri.stim.visContrast>=minMaxVis(1) & blk.tri.stim.visContrast<minMaxVis(2));
end

triType = blk.tri.trialType;
stim = blk.tri.stim;
switch lower(subsetTag)
    case 'blnk'; blk =  prc.filtBlock(blk, triType.blank);
    case 'vl'; blk =  prc.filtBlock(blk, triType.visual & stim.visInitialAzimuth<0);
    case 'vr'; blk =  prc.filtBlock(blk, triType.visual & stim.visInitialAzimuth>0);
    case 'al'; blk =  prc.filtBlock(blk, triType.auditory & stim.audInitialAzimuth<0);
    case 'ar'; blk =  prc.filtBlock(blk, triType.auditory & stim.audInitialAzimuth>0);
    case 'cohl'; blk =  prc.filtBlock(blk, triType.coherent & stim.visInitialAzimuth<0);
    case 'cohr'; blk =  prc.filtBlock(blk, triType.coherent & stim.visInitialAzimuth>0);
    case 'conl'; blk =  prc.filtBlock(blk, triType.conflict & stim.visInitialAzimuth<0);
    case 'conr'; blk =  prc.filtBlock(blk, triType.conflict & stim.visInitialAzimuth>0);
end
end
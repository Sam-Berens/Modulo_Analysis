function [Data,pCover] = getTpatterns_RoiRes(G,subjectId,roiId)
% GETTPATTERNS_ROIRES  Extract Alpha01 T-image patterns at ROI resolution.
%
%   [Data, pCover] = GETTPATTERNS_ROIRES(G, subjectId, roiId)
%
%   Loads a fixed set of first-level SPM T-statistic images for a subject
%   (typically 12 images: two places × six integers) and samples them into
%   the voxel grid of a native-space ROI mask ("ROI resolution"). The
%   output is a matrix of T-values with one row per ROI voxel and one
%   column per T-image/condition.
%
%   Sampling is performed using SPM's volumetric sampling functions, so ROI
%   and T-images do not need to have identical voxel grids, provided they
%   are in the same physical space (i.e., consistent affine transforms).
%
%   INPUTS
%   ------
%   G : char | string
%       Group identifier used by the data pipeline (e.g., 'G0', 'G1').
%
%   subjectId : char | string
%       Subject identifier corresponding to a subject folder in ../../Data.
%
%   roiId : char | string
%       ID of the ROI. The ROI mask is expected to be located via the
%       pipeline helper (e.g., getNativeRoiPath).
%
%   OUTPUTS
%   -------
%   Data : double [nVox × nT]
%       Matrix of T-values sampled at ROI resolution. Rows correspond to
%       ROI voxels (after any optional coverage censoring), columns
%       correspond to T-images/conditions.
%
%   pCover : double
%       Proportion of ROI voxels considered "covered" by the EPI mask.
%
%   NOTES
%   -----
%   - ROI voxels may optionally be censored based on EPI coverage; if so,
%     the returned Data rows correspond to the retained subset of ROI
%     voxels.
%   - Requires SPM on the MATLAB path (spm_vol, spm_read_vols,
%     spm_sample_vol).
%   - Assumes a fixed naming scheme for T-images and ROI files as defined
%     by helper functions in the analysis pipeline.
%
%   EXAMPLE
%   -------
%   [Data, pCover] = getTpatterns_RoiRes('G1','eade18a5','lHippAnt');
%
%   See also GETTIMGS, SPM_VOL, SPM_READ_VOLS, SPM_SAMPLE_VOL.

%% Make sure subjectId is a char
if iscategorical(subjectId)
    subjectId = char(subjectId);
end

%% Get the ROI mask
roiMask.name = getNativeRoiPath(G,subjectId,roiId);
roiMask.V = spm_vol(roiMask.name);
roiMask.M = spm_read_vols(roiMask.V);
roiMask.size = size(roiMask.M);
roiMask.idx = find(roiMask.M > 0.5);

% Convert the idx to mm XYZ:
[x,y,z] = ind2sub(roiMask.size,roiMask.idx);
Xyz = [x,y,z,ones(numel(roiMask.idx),1)]';
Xyz = roiMask.V.mat * Xyz;
roiMask.Xyz = Xyz;
clear Xyz;

%% Get EPI mask
epiMask = getEpiMask(subjectId);

%% Sample the EPI mask for each voxel in the ROI
Vx = epiMask.V.mat \ roiMask.Xyz;
roiCoverage = ...
    spm_sample_vol(epiMask.V,...
    Vx(1,:),Vx(2,:),Vx(3,:),-5)';

%% Get the t-images
Timgs = getTimgs(subjectId); % Note we don't mask here

%% Sample from the t-images
nT = numel(Timgs.P);
Data = nan(numel(roiMask.idx),nT);
for iT = 1:nT
    Vx = Timgs.V(iT).mat \ roiMask.Xyz;
    Data(:,iT) = spm_sample_vol(Timgs.V(iT),Vx(1,:),Vx(2,:),Vx(3,:),-5)';
end

%% Censor Data with roiCoverage and compute pCover
s = roiCoverage>0.5;
Data = Data(s,:);
pCover = mean(s);
return
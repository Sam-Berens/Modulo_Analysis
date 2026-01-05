function [Data,pCover] = getTpatterns_RoiRes(G,subjectId,roiId)

% Data is a [nVox,nPattern] matrix of tStatistics

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
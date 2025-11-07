function [Data] = getNativeData(maskPath,epiPaths,roiPath)

%% Load in volume headers and data
mask.V = spm_vol(maskPath);
epi.V = cellfun(@spm_vol,epiPaths);
roi.V = spm_vol(roiPath);
roi.M = spm_read_vols(roi.V);

%% Find the mm coords (Xyz) for all voxels in ROI
S = roi.M > 0.5;
I = find(S);
[x,y,z] = ind2sub(size(S),I);
Xyz = [x,y,z,ones(size(I))]';
Xyz = roi.V.mat*Xyz;

%% Sample from EPI mask using Xyz
Vx = mask.V.mat\Xyz;
Coverage = ...
    spm_sample_vol(mask.V,...
    Vx(1,:),Vx(2,:),Vx(3,:),-5);

%% Sample from EPI volumes (over time) using Xyz
nTime = numel(epi.V);
Data = nan(nTime,numel(I));
for iTime = 1:nTime
    Vx = epi.V(iTime).mat\Xyz;
    Data(iTime,:) = ...
        spm_sample_vol(epi.V(iTime),...
        Vx(1,:),Vx(2,:),Vx(3,:),-5);
end

%% Censor Data with Coverage
Data = Data(:,Coverage>0.5);
return 
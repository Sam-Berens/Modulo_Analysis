function [Data] = getTpatterns_RoiRes(G,subjectId,roiId,dirs)

%% Load in volume headers and data
mask.V = spm_vol(dirs.maskPath);
mask.M = spm_read_vols(mask.V);
%%Step where you convert roiId into and actual roi
roiPath = getNativeRoiPath(G,subjectId,roiId); 
roi.V = spm_vol(roiPath);

%loop through regressors to build a t-statistic volume which is shaped with
%stimulus in the 4th dim instead of time
nStims = 6;
ts2Collect = cell(nStims,1);
for ii=0:5
    rgName = ['i',num2str(ii)];
    tFp = [dirs.Alpha00 ,filesep,rgName,filesep,'spmT_0001.nii'];
    ts2Collect{(ii+1),1} = tFp;
end
t.V = spm_vol(ts2Collect);

%% Find the mm coords (Xyz) for all voxels in ROI
S = roi.M > 0.5;
I = find(S);
[x,y,z] = ind2sub(size(S),I);
Xyz = [x,y,z,ones(size(I))]';
Xyz = roi.V.mat*Xyz;
%this makes the voxel cordinates where the roi is, in mm space
Vx = t.V{1,1}.mat\Xyz;

%sample the epiMask in roi space to check coverage
Coverage = ...
    spm_sample_vol(mask.V,...
    Vx(1,:),Vx(2,:),Vx(3,:),-5);

%% Sample from t Statistics over regressors using Xyz
nStims = 6;
Data = nan(nStims,numel(I));
for iStim = 1:nStims
    cTs = ...
        spm_sample_vol(t.V{iStim,1},...
        Vx(1,:),Vx(2,:),Vx(3,:),-5);
    %remove any values outside the epiMask
    Data(iStim,:) = cTs(mask.M);
end

%% Censor Data with Coverage
Data = Data(:,Coverage>0.5);
return

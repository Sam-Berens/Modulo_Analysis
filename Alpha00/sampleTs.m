function [Data] = sampleTs(dirs,roiPath)

%% Load in volume headers and data
roi.V = spm_vol(roiPath);
roi.M = spm_read_vols(roi.V);

%loop through regressors to build a t-statistic volume which is shaped with
%stimulus in the 4th dim instead of time
nStims = 6;
ts2Collect = cell(nStims,1);
for ii=0:5
    rgName = ['i',num2str(ii)];
    tFp = [dirs.Alpha00 ,filesep,rgName,filesep,'spmT_0001.nii'];
    ts2Collect{(ii+1),1} = tFp;
end
% ts2Collect = permute(ts2Collect,[2,3,4,1]); %adding singleton dimensions so that we
%eventually end up with one big image matrix where the 4th dim is stimulus
%condition
t.V = spm_vol(ts2Collect);

%% Find the mm coords (Xyz) for all voxels in ROI
S = roi.M > 0.5;
I = find(S);
[x,y,z] = ind2sub(size(S),I);
Xyz = [x,y,z,ones(size(I))]';
Xyz = roi.V.mat*Xyz;

%% Sample from t Statistics over regressors using Xyz
Data = nan(nStims,numel(I));
for iStim = 1:nStims
    Vx = t.V{iStim,1}.mat\Xyz;
    Data(iStim,:) = ...
        spm_sample_vol(t.V{iStim,1},...
        Vx(1,:),Vx(2,:),Vx(3,:),-5);
end

%Had a check to see what it looks like and looks ok but very patchy?
%check if thats normal
% test = nan(size(roi.M));
% test(I) = Data(:,1);
return 

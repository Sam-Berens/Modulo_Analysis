function [Data,roiCoverage] = getTpatterns_EpiRes(G,subjectId,roiId)

dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.Subject = [dirs.Data,filesep,subjectId];
dirs.EPI = [dirs.Subject,filesep,'EPI'];
dirs.Alpha00 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha00'];


% Get ROI mask
roiPath = getNativeRoiPath(G,subjectId,roiId);
roiMask.V = spm_vol(roiPath);
roiMask.M = spm_read_vols(roiMask.V);
roiMask.size = size(roiMask.M);

% Spec the ROI coverage in ROI space
roiCoverage = zeros(size(roiMask.M));

% Get EPI mask
epiMask.name = sprintf('%s%s_%s_epiMask00.nii',dirs.EPI,filesep,subjectId);
epiMask.V = spm_vol(epiMask.name);
epiMask.M = spm_read_vols(epiMask.V);
epiMask.idx = find(epiMask.M > 0.5);
epiMask.size = size(epiMask.M);

% Spec the ROI mask in EPI space!
epiRoi.M = zeros(size(epiMask.M));

% Loop through the epiMask to:
%   1) populate epiRoi by samppling from roiMask
%   2) populate roiCoverage by keeping track of where we have sampled from
for iStim = 1:numel(epiMask.idx)

    % Step 1:
    [vx,vy,vz] = ind2sub(epiMask.size,epiMask.idx(iStim));
    vxEpi = [vx,vy,vz,1]';
    mm = epiMask.V.mat * vxEpi;
    vxRoi = roiMask.V.mat \ mm;
    sample = spm_sample_vol(roiMask.V, ...
        vxRoi(1),vxRoi(2),vxRoi(3),-5);
    epiRoi.M(vx,vy,vz) = sample;

    % Step 2:
    if sample > 0.5
        [w,idxToWeight] = getLerpW(vxRoi(1:3)',roiMask.size);
        roiCoverage(idxToWeight) = roiCoverage(idxToWeight) + w;
        % vxRoi = vxRoi(1:3)';
        % Bvx = getBoundVx(vxRoi);
        % idxToWeight = sub2ind(roiMask.size,Bvx(:,1),Bvx(:,2),Bvx(:,3));
        % dists = pdist2(Bvx,vxRoi);
        % w = (1./dists) ./ sum(1./dists);
        % roiCoverage(idxToWeight) = roiCoverage(idxToWeight) + w;
    end
end

% Get the masked EPI voxel indecies
epiRoi.idx = find(epiRoi.M>0.5);

% Get a matrix of t-stats
Data = nan(numel(epiRoi.idx),6);
for iStim = 1:6
    stimId = ['i',num2str(iStim-1)];
    tFn = [dirs.Alpha00,filesep,stimId,filesep,'spmT_0001.nii'];
    tV = spm_vol(tFn);
    tM = spm_read_vols(tV);
    Data(:,iStim) = tM(epiRoi.idx);
end

heatmap(squareform(pdist(Data','correlation')));

if 0
    %% Load in volume headers and data
    mask.V = spm_vol(dirs.Mask);
    mask.M = spm_read_vols(mask.V);
    roiPath = getNativeRoiPath(G,subjectId,roiId);
    roiMask.V = spm_vol(roiPath);

    % build a 4D vol of t-statistics where 4th dim is stim number
    nStims = 6;
    ts2Collect = cell(nStims,1);
    for iStim=0:5
        rgName = ['i',num2str(iStim)];
        tFp = [dirs.Alpha00 ,filesep,rgName,filesep,'spmT_0001.nii'];
        ts2Collect{(iStim+1),1} = tFp;
    end
    %prepare input to spm_vol such that result (M)
    % is [x, y, z, #Stim]
    ts2Collect = permute(ts2Collect ,[2,3,4,1]);

    %% Find the mm coords (Xyz) for all voxels in EPI
    t.V = spm_vol(ts2Collect);
    t.M = cellfun(@(s)spm_read_vols(s),t.V,'UniformOutput',false);
    t.M = cell2mat(t.M);

    %find voxel coords of t stat M (epi res)
    sz = size(t.M,[1,2,3]); %ignore 4th dim (#Stim) for now
    maxIdx = sz(1) * sz(2) *sz(3);
    [x,y,z] = ind2sub(sz,1:maxIdx);
    % add ones in 4th dim becasue affines are 4d
    Xyz = [x;y;z;ones(1,maxIdx)];
    %transform voxel cords to mm space
    Xyz = t.V{1,1,1,1}.mat*Xyz;


    %% Sample from t Statistics over regressors using Xyz
    nStims = 6;
    Data = cell(nStims,1);

    %Take the xyz coords of voxels(of EPI) ...
    %in mm and convert them into roi voxel space to see if they are in the
    %mask (value is 1) or not(value is 0) at that location
    Vx = roiMask.V.mat\Xyz;
    isEpiInRoi = ...
        spm_sample_vol(roiMask.V,...
        Vx(1,:),Vx(2,:),Vx(3,:),-5);
    isEpiInRoi(isnan(isEpiInRoi)) = 0;
    isEpiInRoi = isEpiInRoi >0.2;%maybe try this as 0.5 looks a bit sparse?

    for iStim = 1:nStims
        %syphon off data for current stimuli
        cTs = t.M(:,:,:,iStim);
        %mask out voxels not in the epi mask
        cTs = cTs.* mask.M;
        %mask out voxels not in the roi mask
        cTs = cTs(isEpiInRoi);
        %remove Ts outside of epiMask coverage
        Data{iStim,1} = cTs;
    end
    Data = cell2mat(Data);
end
return
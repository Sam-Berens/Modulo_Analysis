function [Data,cIsEpiInRoi] = getTpatterns_EpiRes(G,subjectId,roiId,dirs) 

%Data is a [nStim,nVox,nPos] matrix of tStatistics;

%% Load in volume headers and data
mask.V = spm_vol(dirs.maskPath);
mask.M = spm_read_vols(mask.V);
%%Step where you convert roiId into and actual roi
roiPath = getNativeRoiPath(G,subjectId,roiId);
roi.V = spm_vol(roiPath);

%loop through regressors to build a t-statistic volume which is shaped with
%stimulus in the 4th dim instead of time
nStims = 6;
positions = 'ab';
nPos = 2;
for iPos = 1:2
    cPos = positions(iPos);
    ts2Collect = cell(nStims,1);
    for ii=0:5
        rgName = [cPos,num2str(ii)];
        tFp = [dirs.Alpha01 ,filesep,rgName,filesep,'spmT_0001.nii'];
        ts2Collect{(ii+1),1} = tFp;
    end

    %put input to spm_vol into a format such that the resulting...
    % M when its read will be x by y by z by condition
    ts2Collect = permute(ts2Collect ,[2,3,4,1]);

    %% Find the mm coords (Xyz) for all voxels in EPI
    t.V = spm_vol(ts2Collect);
    %this should also have 4th dim being stim, 5th dim being position
    t.M = cellfun(@(s)spm_read_vols(s),t.V,'UniformOutput',false);
    t.M = cell2mat(t.M);

    %% Find the mm coords (Xyz) for all voxels in ROI
    sz = size(t.M,[1,2,3]);
    maxLinInd = sz(1) * sz(2) *sz(3);
    %this gives us all voxel coords of t stat volume (epi res)
    [x,y,z] = ind2sub(sz,1:maxLinInd); %this adds ones in 4th dim becasue affines are 4d
    %this gives us the mm coords of all voxels in the epi image
    Xyz = [x;y;z;ones(1,maxLinInd)];
    Xyz = t.V{1,1,1,1}.mat*Xyz; %doesnt matter which condition we get the .mat from

    %% Sample from t Statistics over regressors using Xyz
    nStims = 6;
    cData = cell(nStims,1);

    %this step says take the xyz coords of voxels(of EPI) ...
    %in mm and put them into voxel 'indices' (in roi voxel space)
    %rember they aren't whole indices hence why we're interpolating
    Vx = roi.V.mat\Xyz;
    %now we want to check to see is there a 1 in the roi mask for each
    %point of epi space
    cIsEpiInRoi = ...
        spm_sample_vol(roi.V,...
        Vx(1,:),Vx(2,:),Vx(3,:),-5);
    cIsEpiInRoi(isnan(cIsEpiInRoi)) = 0;
    cIsEpiInRoi = cIsEpiInRoi >0.5;

    for iStim = 1:nStims
        %syphon off current data
        cTs = t.M(:,:,:,iStim);
        %mask out bits not in the epi mask;
        cTs = cTs.* mask.M;
        %...then we pick out the t data at the point where there are ones...
        cTs = cTs(cIsEpiInRoi); %here we are linearly indexing
        %remove Ts outside of epiMask coverage
        cData{iStim,1} = cTs; %here we are lin indexing using ...
        % where the mask is nonzero (should only be ones or zeros)...
        %also we're doing it like this to flatten the data
    end

    if strcmp(cPos,'a')
        %if its the first position then build Data for first time
        % and put current data in the first layer of the 3rd Dim
        nVox = size(cData{1,1},2);
        Data = nan(nStims,nVox,nPos);
        Data(:,:,1) = cell2mat(cData);
        isEpiInRoi(:,1) = cIsEpiInRoi;
    else
        %if its pos b we're filling in the 2nd layer of the 3rd Dim
        Data(:,:,2) = cell2mat(cData); 
        isEpiInRoi(:,2) = cIsEpiInRoi;
    end
end
return
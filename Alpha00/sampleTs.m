function [Data] = sampleTs(dirs,roiPath,maskPath,m)

%% Load in volume headers and data
mask.V = spm_vol(maskPath);
mask.M = spm_read_vols(mask.V);
roi.V = spm_vol(roiPath);
roi.M = spm_read_vols(roi.V);

%loop through regressors to build a t-statistic volume which is shaped with
%stimulus in the 4th dim instead of time
nStims = 6;
nPos = 2;
ts2Collect = cell(nStims,nPos);
%%COMENBACK AND PUT IN NEW DIM for position
for ii=0:5
    rgName = ['i',num2str(ii)];
    tFp = [dirs.Alpha00 ,filesep,rgName,filesep,'spmT_0001.nii'];
    ts2Collect{(ii+1),1} = tFp;
end

t.V = spm_vol(ts2Collect);

%m = the interpolation method you want to use for sampling
if strcmp(m,'push')
    [Data] = interp1(roi,t,mask);
elseif strcmp(m,'pull')
    [Data] = interp2(ts2Collect,t,mask);
end

%Had a check to see what it looks like and looks ok but very patchy?
%check if thats normal
% test = nan(size(roi.M));
% test(I) = Data(:,1);
return 




function [Data] = interp1(roi,t,mask)
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
    cTs = ...
        spm_sample_vol(t.V{iStim,1},...
        Vx(1,:),Vx(2,:),Vx(3,:),-5);
    %remove any values outside the epiMask
    Data(iStim,:) = cTs(mask.M);
end
return

function [Data] = interp2(ts2Collect,t,mask)
%% Find the mm coords (Xyz) for all voxels in EPI
% S = roi.M > 0.5;
% I = find(S);
%put V into a format such that the resulting M will be x by y by z by
%condition
ts2Collect = permute(ts2Collect ,[2,3,4,1]);
t.V = spm_vol(ts2Collect);
%this should also have 4th dim being condition
t.M = spm_read_vols(t.V);
%this just gives us all voxel coords of t stat volume (epi res)
[x,y,z] = ind2sub(size(t.M,[1,2,3]));
Xyz = [x,y,z,ones(size(t.M,[1,2,3]))]'; %this adds ones in 4th dim becasue affines are 4d
%this gives us the mm coords of all voxels in the epi image
Xyz = t.V(1).mat*Xyz;

%% Sample from t Statistics over regressors using Xyz
Data = nan(nStims,numel(I));
for iStim = 1:nStims
    %this step says take the xyz coords of voxels(of EPI) ...
    %in mm and put them into voxel 'indices' (in roi voxel space)
    %rember they aren't whole indices hence why we're interpolating
    Vx = roi.V.mat\Xyz;
    %now we want to check to see is there a 1 in the roi mask for each
    %point of epi space
    isEpiInRoi = ...
        spm_sample_vol(roi.V,...
        Vx(1,:),Vx(2,:),Vx(3,:),0); %check which type we want to do... this is nearest neighbour
    %syphon off current data
    cTs = t.M(:,:,:,iStim);
     %...then we pick out the t data at the point where there are ones...
    cTs = cTs(isEpiInRoi); %here we are linearly indexing
    %remove Ts outside of epiMask coverage
    Data(iStim,:) = cTs(find(mask.M)); %here we are lin indexing using ...
    % where the mask is nonzero (should only be ones or zeros)...
    %also we're doing it like this to flatten the data
end
return
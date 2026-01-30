function [] = zX6_makeEpsilon()
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
subjectIds = char(subjectIds);
for iSubject=1:nSubs
    subjectId = subjectIds(iSubject,:);
    makeEpsilon(subjectId,G);
end
return


function [] = makeEpsilon(subjectId,G)
dirs.Data = '../../Data';
dirs.Subject = [dirs.Data,filesep,char(subjectId)];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
dirs.Mdl03 = [dirs.Alpha01,filesep,'Mdl03']; %TO DO DECIDE FOLDER NAME
dirs.G1 = [dirs.Mdl03,filesep,G];
if ~exist(dirs.Mdl03 ,'dir')
    mkdir(dirs.Mdl03);
end

epiMask = getEpiMask(subjectId);

%% Get a matrix of t-stats for all conditions
%remember that a stims are stacked ontop of b stim
[tImgs] = getTimgs(subjectId,epiMask);

yDepth = 8; %Y is going to have 8 outputs (1st 4 for coloc=-1 2nd 4 for coloc=+1)
%% loop through searchlight centres to test produce epsilon term from precursor to mdl03
r = 3; 
Y = searchlight3D(3,@qFunc_Mdl2,epiMask,tImgs,yDepth); %reminder N is the volume of each searchlight

%we're going to save all stats maps as nifti so we can norm q images
%to mni space and so that we can do QA checks on the others
coLocM1 = Y(:,:,:,1);
coLocP1 = Y(:,:,:,2);

coLocM1Mask = ~isnan(coLocM1);
coLocP1Mask = ~isnan(coLocP1);


%for the q analysis we set nans to be 1s? so that when they got
%interpolated we were only biasing the images towards q = 1? because nans
%get turned into zeros in the interpolation?

images = {coLocM1,coLocP1,coLocM1Mask,coLocP1Mask};
names = {'coLocM1','coLocP1','coLocM1Mask','coLocP1Mask'};
descrip = sprintf('Epsilon map for nonlin RDM model, searchlight r=%i',r);
%save to their mdl03 folder
for iIm=1:numel(images)
    im.M = images{iIm};
    im.V = epiMask.V;
    im.V.dt(1) = 64;
    label = names{iIm};
    im.V.fname = [dirs.Mdl03,filesep,label,'.nii'];
    im.V.descrip = descrip;
    for iSlice=1:2
        V = im.V;
        V.n = [iSlice,1];
        M = im.M(:,:,:,iSlice);
        % this looks like it overwrites it but
        % it just adds the slice to any existing file
        spm_write_vol(V,M); 
    end
end
return
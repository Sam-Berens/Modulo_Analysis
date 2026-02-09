function [] = zX03_makeQ()
% Cd out
wd = pwd;
cd ..;

G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
subjectIds = char(subjectIds);
for iSubject=1:nSubs
    subjectId = subjectIds(iSubject,:);
    makeQ(subjectId,G);
end
% Cd back
cd(wd);
return


function [] = makeQ(subjectId,G)
dirs.Data = '../../Data';
dirs.Subject = [dirs.Data,filesep,char(subjectId)];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
dirs.Mdl01 = [dirs.Alpha01,filesep,'Mdl01'];
dirs.G1 = [dirs.Mdl01,filesep,G];
if ~exist(dirs.Mdl01 ,'dir')
    mkdir(dirs.Mdl01);
end

epiMask = getEpiMask(subjectId);

%% Get a matrix of t-stats for all conditions
%remember that a stims are stacked ontop of b stim
[tImgs] = getTimgs(subjectId,epiMask);

yDepth = 4; %this is because the axis model has 3 params plus the error term
%% loop through searchlight centres to test produce q term from precursor to mdl01
r = 3; 
Y = searchlight3D(3,@qFunc,epiMask,tImgs,yDepth); %reminder N is the volume of each searchlight

%we're going to save all stats maps as nifti so we can norm q images
%to mni space and so that we can do QA checks on the others
B0 = Y(:,:,:,1);
B1 = Y(:,:,:,2);
Q = Y(:,:,:,3);
%error of the model
err = Y(:,:,:,4);
%images needed for model 1
lQ = log(Q);
qMask = ~isnan(Q);

images = {Q,B0,B1,err,lQ,qMask};
names = {'q','b0','b1','error','lQ','qMask'};
descrip = sprintf('Stat map for nonlin RDM model, searchlight r=%i',r);
for iIm=1:numel(images)
    im.M = images{iIm};
    im.V = epiMask.V;
    im.V.dt(1) = 64;
    label = names{iIm};
    %save to their mdl01 folder
    im.V.fname = [dirs.Mdl01,filesep,label,'.nii'];
    im.V.descrip = descrip;
    spm_write_vol(im.V,im.M);
end
return

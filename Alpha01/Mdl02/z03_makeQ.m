function [] = z03_makeQ()
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
% Cd back in
cd(wd);
return


function [] = makeQ(subjectId,G)
dirs.Data = '../../Data';
dirs.Subject = [dirs.Data,filesep,char(subjectId)];
dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
dirs.Mdl02 = [dirs.Alpha01,filesep,'Mdl02'];
dirs.G1 = [dirs.Mdl02,filesep,G];
if ~exist(dirs.Mdl02 ,'dir')
    mkdir(dirs.Mdl02);
end

epiMask = getEpiMask(subjectId);

%% Get a matrix of t-stats for all conditions
%remember that a stims are stacked ontop of b stim
[tImgs] = getTimgs(subjectId,epiMask);

yDepth = 8; % Y is going to have 8 outputs (1st 4 for coloc=-1 2nd 4 for coloc=+1)
%% loop through searchlight centres to test produce q term from precursor to mdl01
r = 3;
Y = searchlight3D(r,@qFunc2,epiMask,tImgs,yDepth); %reminder N is the volume of each searchlight
B0 = cat(4,Y(:,:,:,1),Y(:,:,:,5));
B1 = cat(4,Y(:,:,:,2),Y(:,:,:,6));
Q = cat(4,Y(:,:,:,3),Y(:,:,:,7));
%error of the model
err = cat(4,Y(:,:,:,4),Y(:,:,:,8));
%images needed for model 2
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
    im.V.fname = [dirs.Mdl02,filesep,label,'.nii'];
    im.V.descrip = descrip;
    for iSlice=1:2
        V = im.V;
        V.n = [iSlice,1];
        M = im.M(:,:,:,iSlice);
        % this looks like it overwrites the .nii but
        % it just adds the slice to any existing file
        spm_write_vol(V,M); 
    end
end

return

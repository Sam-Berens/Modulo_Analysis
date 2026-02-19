function [] = z03_runSearchlight(G)
wd = pwd;
%Back up into Alpha01
cd ..
% Get the list of subjectIds
subjectIds = getSubjectIds(G);

% Set the path to the /Data directory
dirs.Data = ['..',filesep,'..',filesep,'Data'];

% Loop through subjects
for iSubject = 1:numel(subjectIds)

    % Set the current subject ID
    cSid = subjectIds(iSubject);

    % Get the EPI mask (native space)
    Mask = getEpiMask(cSid);

    % Get the Timags (native space)
    Timgs = getTimgs(cSid);

    % Set the output directory
    dirs.Subject = [dirs.Data,filesep,char(cSid)];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
    dirs.Mdl05 = [dirs.Alpha01,filesep,'Mdl05'];
    if ~exist(dirs.Mdl05,'dir')
        mkdir(dirs.Mdl05);
    end

    % Set the colocation = -1 output header
    Vn = Mask.V;
    Vn.fname = [dirs.Searchlight,filesep,'zTemplate_colocation=-1.nii'];
    Vn.dt(1) = 64;
    Vn.descrip = 'Alpha01 zTemplate searchlight: colocation=-1;';

    % Set the colocation = +1 output header
    Vp = Mask.V;
    Vp.fname = [dirs.Searchlight,filesep,'zTemplate_colocation=+1.nii'];
    Vp.dt(1) = 64;
    Vp.descrip = 'Alpha01 zTemplate searchlight: colocation=+1;';

    % Set the visSim output header
    Vv = Mask.V;
    Vp.fname = [dirs.Searchlight,filesep,'zTemplate_visSim.nii'];
    Vp.dt(1) = 64;
    Vp.descrip = 'Alpha01 zTemplate searchlight: visSim;';

    % Set the neighborhood size image
    Vc = Mask.V;
    Vc.fname = [dirs.Searchlight,filesep,'hoodSize.nii'];
    Vc.dt(1) = 16;
    Vc.descrip = 'Alpha01 zTemplate searchlight: neighborhood size (#vx);';

    %this is a beta for coloc-1, coloc+1 and VisStruct
    Ydepth = 3;
    %TO DO _ PICK RAD
    r = 5;
    % Run the searchlight
    [Z,N] = searchlight3D(...
        r,...         Radius   -
        @searchFun,... Function
        Mask,...       Mask
        Timgs,...      Data
        Ydepth,...     Output depth
        char(cSid));

    % Save the results
    Mn = Z(:,:,:,1);
    Mp = Z(:,:,:,2);
    Mv = Z(:,:,:,3);
    spm_write_vol(Vn,Mn);
    spm_write_vol(Vp,Mp);
    spm_write_vol(Vv,Mv);
    spm_write_vol(Vc,N);
end
%return to Mdl04 folder
cd(wd);
return


function [z] = searchFun(M)

persistent H_ h S;
if isempty(S) || ...
        isempty(h) || ...
        isempty(H_)
    % Construct the basic 6x6 similarity hypothesis
    simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
    H_ = nan(6);
    for ii = 1:36
        [x,y] = ind2sub([6,6],ii);
        H_(ii) = simFun(x,y);
    end
    h = H_(tril(true(6),-1));

    %make the selectors for the real data RSM
    S.n.ab = kron([0,0;1,0],tril(true(6), -1));
    S.n.ba = kron([0,0;1,0],trih(true(6), 1));
    S.p.a = kron([1,0;0,0],tril(true(6), -1));
    S.p.b =  kron([0,0;0,1],tril(true(6), -1));
end
%get the 12x12 RSM for the searchlight image
R = corr(M);

%load in predicted visual similarity based on their permutation
%and layer X of denseNet-169
%TO DO - just run this once and save ! and then load
V = getVisCorr();

%% Permute the values of V to reflect the images assigned to each number for
%% that subject
%this assumes that perm works such that e.g. image i03 becomes
%image 1, i.e. the index is indicated by the position, not the number
imgPerm = getImagePerm(subjectId);
%get rid of zero ordering (currently rsm is arranged such that column 1 is
%comparisons with image i00.png, so we want the 0 in imgPerm to get mapped
%to column 1 (before all the permuting)
imgPerm = imgPerm + 1;
idxs = (1:36)';
[x,y] = ind2sub(size(V),idxs);
[~,pidx] = sort(imgPerm);
nx = arrayfun(@(x) pidx(x),x);
ny = arrayfun(@(y) pidx(y),y);
prmV = nan(6);
for ii=1:36
    prmV(nx(ii),ny(ii)) = V(idxs(ii));
end
%select out lower tri
v = prmV(tril(true(6),-1));

%scale all vectors such that they are projecting from the origin
%and are zscored
h = zscore(h);
h = repmat(h,[2,1]);
v = zscore(v);
v = repmat(v,[4,1]);

%% prep y column of design mat
%select out colcoation based vals of R and zscore within aa/bb/ba comp type
yn = zscore([R(S.n.ab);R(S.n.ba)]);
yp.a = zscore(R(S.p.a));
yp.b = zscore(R(S.p.b));
y = [yn;yp.a;yp.b];

%construct the design matrix
X = [kron(eye(2),h),v];
%do regression and Fischer-transform the results
z = atanh(pinv(X)*y);
return
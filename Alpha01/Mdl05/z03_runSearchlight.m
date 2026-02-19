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

    %Get subject's image perm, undo zero-ordering, get pidx for vis RSM
    imgPerm = getImagePerm(char(cSid));
    imgPerm = imgPerm + 1;
    [~,pidx] = sort(imgPerm);

    % Set the output directory
    dirs.Subject = [dirs.Data,filesep,char(cSid)];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
    dirs.Mdl05 = [dirs.Alpha01,filesep,'Mdl05'];
    if ~exist(dirs.Mdl05,'dir')
        mkdir(dirs.Mdl05);
    end

    %TO DO _ PICK RAD
    r = 5;
    % Set the colocation = -1 output header
    Vn = Mask.V;
    Vn.fname = [dirs.Mdl05,filesep,'zTemplate_colocation=-1.nii'];
    Vn.dt(1) = 64;
    Vn.descrip = sprintf('Alpha01 zTemplate searchlight: colocation=-1, r=%i',r);

    % Set the colocation = +1 output header
    Vp = Mask.V;
    Vp.fname = [dirs.Mdl05,filesep,'zTemplate_colocation=+1.nii'];
    Vp.dt(1) = 64;
    Vp.descrip = sprintf('Alpha01 zTemplate searchlight: colocation=+1, r=%i',r);

    % Set the visSim output header
    Vv = Mask.V;
    Vp.fname = [dirs.Mdl05,filesep,'zTemplate_visSim.nii'];
    Vp.dt(1) = 64;
    Vp.descrip = sprintf('Alpha01 zTemplate searchlight: visSim, r=%i',r);

    % Set the neighborhood size image
    Vc = Mask.V;
    Vc.fname = [dirs.Mdl05,filesep,'hoodSize.nii'];
    Vc.dt(1) = 16;
    Vc.descrip = sprintf(...
        'Alpha01 zTemplate searchlight: neighborhood size (#vx), r=%i',r);

    %this is a beta for coloc-1, coloc+1 and VisStruct
    Ydepth = 3;
    % Run the searchlight
    [Z,N] = searchlight3D(...
        r,...         Radius   -
        @(M) searchFun(M,pidx),... Function
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


function [z] = searchFun(M,pidx)

persistent H_ h S v;
if isempty(v) ||...
        isempty(S) || ...
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
    S.n.ab = logical(kron([0,0;1,0],tril(true(6),-1)));
    S.n.ba = logical(kron([0,0;1,0],triu(true(6),1)));
    S.p.a = logical(kron([1,0;0,0],tril(true(6),-1)));
    S.p.b = logical(kron([0,0;0,1],tril(true(6),-1)));

    %load in predicted visual similarity based on their permutation
    %and layer X of denseNet-169
    V = getVisCorr();
    %% Permute the values of V to reflect the images assigned to each number for
    %% that subject
    %this assumes that perm works such that e.g. image i03 becomes
    %image 1, i.e. the index is indicated by the position, not the number
    %get rid of zero ordering (currently rsm is arranged such that column 1 is
    %comparisons with image i00.png, so we want the 0 in imgPerm to get mapped
    %to column 1 (before all the permuting)
    idxs = (1:36)';
    [x,y] = ind2sub(size(V),idxs);
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

end
%get the 12x12 RSM for the searchlight image
R = corr(M);

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
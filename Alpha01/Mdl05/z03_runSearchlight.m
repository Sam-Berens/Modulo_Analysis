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

    %Get subject's image perm, undo zero-ordering
    imgPerm = getImagePerm(char(cSid));
    imgPerm = imgPerm + 1;
    
    % Set the output directory
    dirs.Subject = [dirs.Data,filesep,char(cSid)];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
    dirs.Mdl05 = [dirs.Alpha01,filesep,'Mdl05'];
    if ~exist(dirs.Mdl05,'dir')
        mkdir(dirs.Mdl05);
    end

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
    Vv.fname = [dirs.Mdl05,filesep,'zTemplate_visSim.nii'];
    Vv.dt(1) = 64;
    Vv.descrip = sprintf('Alpha01 zTemplate searchlight: visSim, r=%i',r);

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
       @(M) mdl5Func(M,imgPerm),... Function
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


function [z] = mdl5Func(M,imgPerm)

persistent prev_imgPerm S IX;
if isempty(prev_imgPerm) || isempty(S) || isempty(IX) || ...
    ~isequal(prev_imgPerm, imgPerm)
        
    % Construct the basic 6x6 similarity hypothesis
    simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
    H_ = nan(6);
    for ii = 1:36
        [x,y] = ind2sub([6,6],ii);
        H_(ii) = simFun(x,y);
    end
    h = H_(tril(true(6),-1));

    % Make the selectors for the real data RSM
    S.n.ab = logical(kron([0,0;1,0],tril(true(6),-1)));
    S.n.ba = logical(kron([0,0;1,0],triu(true(6),1)));
    S.p.a = logical(kron([1,0;0,0],tril(true(6),-1)));
    S.p.b = logical(kron([0,0;0,1],tril(true(6),-1)));

    % Load in predicted visual similarity from denseNet-169
    V = getVisCorr();
    
    % Permute the rows and cols of V
    PV = V(imgPerm,imgPerm);
    %select out lower tri
    v = PV(tril(true(6),-1));
    %scale all vectors such that they are projecting from the origin
    %and are zscored
    h = zscore(h);
    h = repmat(h,[2,1]);
    v = zscore(v);
    v = repmat(v,[4,1]);

    % Construct the design matrix (X)
    X = [kron(eye(2),h),v];

    % Invert for regression
    IX = pinv(X);

    % Set/update prev_imgPerm
    prev_imgPerm = imgPerm;

end
%get the 12x12 RSM for the searchlight image
R = corr(M);

%% prep y column of design mat
%select out colcoation based vals of R and zscore within aa/bb/ba comp type
yn = zscore([R(S.n.ab);R(S.n.ba)]);
yp.a = zscore(R(S.p.a));
yp.b = zscore(R(S.p.b));
y = [yn;yp.a;yp.b];

%do regression and Fischer-transform the results
z = atanh(IX*y);
return
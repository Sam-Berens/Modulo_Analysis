function [] = z03_runSearchlight(G)

% Searchlight radius
r = 5;

% Get the list of subjectIds
subjectIds = getSubjectIds(G);

% Move up into Alpha01
wd = pwd;
cd ..;

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

    % Get subject's image perm, undo zero-ordering
    imgPerm = getImgPerm(char(cSid));
    imgPerm = imgPerm + 1;

    % Set the output directory
    dirs.Subject = [dirs.Data,filesep,char(cSid)];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];
    dirs.Mdl05 = [dirs.Alpha01,filesep,'Mdl05'];
    if ~exist(dirs.Mdl05,'dir')
        mkdir(dirs.Mdl05);
    end

    % Set the colocation = -1 output header
    Vn = Mask.V;
    Vn.fname = [dirs.Mdl05,filesep,'zTemplate_colocation=-1.nii'];
    Vn.dt(1) = 64;
    Vn.descrip = sprintf(...
        'Alpha01 Mdl05 searchlight: colocation=-1, r=%i',r);

    % Set the colocation = +1 output header
    Vp = Mask.V;
    Vp.fname = [dirs.Mdl05,filesep,'zTemplate_colocation=+1.nii'];
    Vp.dt(1) = 64;
    Vp.descrip = sprintf(...
        'Alpha01 Mdl05 searchlight: colocation=+1, r=%i',r);

    % Set the visSim output header
    Vv = Mask.V;
    Vv.fname = [dirs.Mdl05,filesep,'zTemplate_visSim.nii'];
    Vv.dt(1) = 64;
    Vv.descrip = sprintf(...
        'Alpha01 Mdl05 searchlight: visSim, r=%i',r);

    % Set the neighborhood size image
    Vc = Mask.V;
    Vc.fname = [dirs.Mdl05,filesep,'hoodSize.nii'];
    Vc.dt(1) = 16;
    Vc.descrip = sprintf(...
        'Alpha01 Mdl05 searchlight: neighborhood size (#vx), r=%i',r);

    % Run the searchlight
    [Z,N] = searchlight3D(...
        r,...     Radius
        @(M)      mdl5Func(M,imgPerm),...
        Mask,...  Mask
        Timgs,... Data
        3,...     Output depth
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

% Return to the working directory
cd(wd);
return

function [z] = mdl5Func(M,imgPerm)
persistent nanzscore VisSim prev_imgPerm S IX;
if isempty(nanzscore) || isempty(VisSim)
    nanzscore = @(v) (v-mean(v,'omitmissing')) ./ std(v,'omitmissing');

    temp = load('VisualSim.mat');
    VisSim = temp.VisualSim;
end

if isempty(prev_imgPerm) || isempty(S) || isempty(IX) || ...
        ~isequal(prev_imgPerm, imgPerm)

    % Construct the basic 6x6 similarity hypothesis
    simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
    H_ = nan(6);
    for ii = 1:36
        [x,y] = ind2sub([6,6],ii);
        H_(ii) = simFun(x,y);
    end

    % Hypothesis for colocation=-1
    Hn = kron(ones(2),H_);
    Hn(logical(kron([1,1;0,1],ones(6)))) = NaN;
    Hn(logical(kron([0,0;1,0],eye(6)))) = NaN;

    Haa = kron(ones(2),H_);
    Haa(logical(kron([0,1;1,1],ones(6)))) = NaN;
    Haa(triu(true(12))) = NaN;

    Hbb = kron(ones(2),H_);
    Hbb(logical(kron([1,1;1,0],ones(6)))) = NaN;
    Hbb(triu(true(12))) = NaN;

    Hp = sum(cat(3,Haa,Hbb),3,'omitmissing');
    Hp(isnan(Haa) & isnan(Hbb)) = NaN;

    % Make the selectors for the real data RSM
    S.n = ~isnan(Hn);
    S.aa = ~isnan(Haa);
    S.bb = ~isnan(Hbb);
    S.p = ~isnan(Hp);
    S.all = S.n | S.p;

    % Load in predicted visual similarity from denseNet-169
    V = VisSim(imgPerm,imgPerm);
    V = kron(ones(2),V);

    % Construct the design matrix
    X = [...
        nanzscore(Hn(S.all)),...
        nanzscore(Hp(S.all)),...
        nanzscore(V(S.all))];

    % Replace N/A values with mean
    X(isnan(X)) = 0;

    % Invert for regression
    IX = pinv(X);

    % Set/update prev_imgPerm
    prev_imgPerm = imgPerm;

end

% Get the 12x12 RSM for the searchlight image
R = corr(M);

% Zscore within aa/bb/ba pairs
R(S.n) = nanzscore(R(S.n));
R(S.aa) = nanzscore(R(S.aa));
R(S.bb) = nanzscore(R(S.bb));
y = R(S.all);

% Do regression and Fischer-transform the results
z = atanh(IX*y);
return
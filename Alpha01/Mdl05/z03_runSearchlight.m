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
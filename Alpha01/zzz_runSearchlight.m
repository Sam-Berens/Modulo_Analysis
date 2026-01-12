function [] = zzz_runSearchlight(G)

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
    dirs.Searchlight = [dirs.Alpha01,filesep,'Searchlight'];
    if ~exist(dirs.Searchlight,'dir')
        mkdir(dirs.Searchlight);
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

    % Set the neighborhood size image
    Vc = Mask.V;
    Vc.fname = [dirs.Searchlight,filesep,'hoodSize.nii'];
    Vc.dt(1) = 2;
    Vc.descrip = 'Alpha01 zTemplate searchlight: neighborhood size (#vx);';

    % Run the searchlight
    [Z,N] = searchlight3D(...
        3,...          Radius
        @searchFun,... Function
        Mask,...       Mask
        Timgs,...      Data
        2,...          Output depth
        char(cSid));

    % Save the results
    Mn = Z(:,:,:,1);
    Mp = Z(:,:,:,2);
    spm_write_vol(Vn,Mn);
    spm_write_vol(Vp,Mp);
    spm_write_vol(Vc,N);
end
return

function [z] = searchFun(X)

persistent Sn Spa Spb hn hpa hpb;
if isempty(Sn) || isempty(Spa) || isempty(Spb) || isempty(hn) || isempty(hpa) || isempty(hpb)
    
    % Construct the basic 6x6 similarity hypothesis
    simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
    H_ = nan(6);
    for ii = 1:36
        [x,y] = ind2sub([6,6],ii);
        H_(ii) = simFun(x,y);
    end

    % Kron to get a 12x12
    H_ = kron(ones(2),H_);

    % Fork colocation = -1 and colocation = +1
    Hn_ = H_;
    Hn_(logical(kron(eye(2),ones(6)))) = NaN;
    Hpa_ = H_;
    Hpa_(logical(kron([0,1;1,1],ones(6)))) = NaN;
    Hpb_ = H_;
    Hpb_(logical(kron([1,1;1,0],ones(6)))) = NaN;

    % Remove the visual effect from Hn
    Hn_(logical(kron([0,1;1,0],eye(6)))) = NaN;

    % Remove upper triangle of H*, including the diagonal
    Hn_(triu(true(12))) = NaN;
    Hpa_(triu(true(12))) = NaN;
    Hpb_(triu(true(12))) = NaN;

    % Set the ~isnan selector
    Sn = ~isnan(Hn_);
    Spa = ~isnan(Hpa_);
    Spb = ~isnan(Hpb_);

    % Vectorise H_ excluding nans
    hn_ = Hn_(Sn);
    hpa_ = Hpa_(Spa);
    hpb_ = Hpb_(Spb);

    % Prep for dot-projuct
    hn_ = hn_ - mean(hn_);
    hn_ = hn_ ./ norm(hn_);
    hn = hn_';

    hpa_ = hpa_ - mean(hpa_);
    hpa_ = hpa_ ./ norm(hpa_);
    hpa = hpa_';

    hpb_ = hpb_ - mean(hpb_);
    hpb_ = hpb_ ./ norm(hpb_);
    hpb = hpb_';
end

R = corr(X);

rn = R(Sn);
rn = rn - mean(rn);
rn = rn ./ norm(rn);
zn = atanh(hn*rn);

rpa = R(Spa);
rpa = rpa - mean(rpa);
rpa = rpa ./ norm(rpa);
zpa = atanh(hpa*rpa);

rpb = R(Spb);
rpb = rpb - mean(rpb);
rpb = rpb ./ norm(rpb);
zpb = atanh(hpb*rpb);

z = [zn;
    (zpa+zpb)/2];
return
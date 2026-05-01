function [] = z02_makeL1PPIRegressors(G,roiInfo)
%% Notes

% roiInfo.fn = '/mnt/Erebus/Modulo/Data/_Group/MniRois/_Cluster-Alpha01-Mdl05a_+zPnonc_rVisual.nii';
% roiInfo.id = 'rVisual';

% roiInfo.fn = '/mnt/Erebus/Modulo/Data/_Group/MniRois/_Cluster-Alpha01-Mdl05a_+zPnonc_lVisual.nii';
% roiInfo.id = 'lVisual';

%ideally we'd want to just construct the fn like this but this doesnt quite
%work for the clusters because we're not using their fullname in later
%labelling:     roiFn = ['/mnt/Erebus/Modulo/Data/_Group/MniRois/',roiId];

%% Takes group and a struct with the fields .fn and .Id which contain the filename and short identifier for the roi
roiFn = roiInfo.fn;
roiId = roiInfo.id;

subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
for iSubject = 1:nSubs
    subjectId = char(subjectIds(iSubject));

    % get subject-specific filenames and info
    dirs.Subject = fullfile('/mnt','Erebus','Modulo','Data',subjectId);
    dirs.A02 = fullfile(dirs.Subject,'Analysis','A02');
    dirs.PPI = fullfile(dirs.A02,sprintf('PPI_%s',roiId));
    if ~exist(dirs.PPI,'dir')
        mkdir(dirs.PPI)
    end
    epiMask = fullfile(dirs.Subject,'EPI',G,sprintf('w_%s_epiMask00.nii',subjectId));
    spmMatFn = fullfile(dirs.A02,'SPM.mat');
    tmp = load(spmMatFn);
    nRuns = size(tmp.SPM.Sess,2);
    clear tmp;

    spmBatch = cell(1,nRuns*3+1);
    jobNumber = 1;
    spmBatch{jobNumber}.spm.stats.con.spmmat = {spmMatFn};
    spmBatch{jobNumber}.spm.stats.con.consess{1}.fcon.name = 'Eye(2)';
    spmBatch{jobNumber}.spm.stats.con.consess{1}.fcon.weights = eye(2);
    spmBatch{jobNumber}.spm.stats.con.consess{1}.fcon.sessrep = 'replsc';
    spmBatch{jobNumber}.spm.stats.con.delete = 1;
    
    for iRun = 1:nRuns

        % Create the VOI
        jobNumber = jobNumber + 1;
        spmBatch{jobNumber}.spm.util.voi.spmmat = {spmMatFn};
        spmBatch{jobNumber}.spm.util.voi.adjust = 1;
        spmBatch{jobNumber}.spm.util.voi.session = iRun;
        spmBatch{jobNumber}.spm.util.voi.name = roiId;
        spmBatch{jobNumber}.spm.util.voi.roi{1}.mask.image = {roiFn};
        spmBatch{jobNumber}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        spmBatch{jobNumber}.spm.util.voi.roi{2}.mask.image = {epiMask};
        spmBatch{jobNumber}.spm.util.voi.roi{2}.mask.threshold = 0.5;
        spmBatch{jobNumber}.spm.util.voi.expression = 'i1&i2';

        %% Set some useful vars
        voiFn = fullfile(dirs.A02,sprintf('VOI_%s_%i.mat',roiId,iRun));
        ppiName.a = sprintf('R%i-%s*a',iRun,roiId);
        ppiName.b = sprintf('R%i-%s*b',iRun,roiId);

        %% Create the PPI regressor for a
        jobNumber = jobNumber + 1;
        spmBatch{jobNumber}.spm.stats.ppi.spmmat = {spmMatFn};
        spmBatch{jobNumber}.spm.stats.ppi.type.ppi.voi = {voiFn};
        % The first element of ppi.u (below) indexes the regressor to use
        spmBatch{jobNumber}.spm.stats.ppi.type.ppi.u = [1 1 1];
        spmBatch{jobNumber}.spm.stats.ppi.name = ppiName.a ;
        spmBatch{jobNumber}.spm.stats.ppi.disp = 0;

        %% Create the PPI regressor for b
        jobNumber = jobNumber + 1;
        spmBatch{jobNumber}.spm.stats.ppi.spmmat = {spmMatFn};
        spmBatch{jobNumber}.spm.stats.ppi.type.ppi.voi = {voiFn};
        % The first element of ppi.u (below) indexes the regressor to use
        spmBatch{jobNumber}.spm.stats.ppi.type.ppi.u = [2 1 1];
        spmBatch{jobNumber}.spm.stats.ppi.name = ppiName.b ;
        spmBatch{jobNumber}.spm.stats.ppi.disp = 0;
    end

    spm_jobman('initcfg');
    spm_jobman('run',spmBatch);
    toMove = [dirs.A02,filesep,'PPI_*-',roiId,'*'];
    movefile(toMove,dirs.PPI);
end
return
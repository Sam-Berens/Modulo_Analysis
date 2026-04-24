
% roiInfo.fn = '/mnt/Erebus/Modulo/Data/_Group/MniRois/_Cluster-Alpha01-Mdl05a_+zPnonc_rVisual.nii';
% roiInfo.id = 'rVisual';

%ideally we'd want to just contruct the fn like this but this doesnt quite
%work for the clusters because we're not using their fullname in later
%labelling:     roiFn = ['/mnt/Erebus/Modulo/Data/_Group/MniRois/',roiId];
function [] = z02_makePPI(G,roiInfo)
%% Takes group and a struct with the fields .fn and .Id which contain the filename and short identifier for the roi
roiFn = roiInfo.fn;
roiId = roiInfo.id;

subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
for iSubject = 1:nSubs
    subjectId = char(subjectIds(iSubject));
    %get subject-specific filenames and info
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

    spmBatch = cell(1,(nRuns*3)+1);
    spmBatch{1}.spm.stats.con.spmmat = {spmMatFn};
    spmBatch{1}.spm.stats.con.consess{1}.fcon.name = 'Eye(2)';
    spmBatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 0
        0 1];
    spmBatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'replsc';
    spmBatch{1}.spm.stats.con.delete = 1;

    batchIdxs = 2:(nRuns*3)+1;
    for iRun = 1:nRuns
        voiFn = fullfile(dirs.A02,sprintf('VOI_%s_%i.mat',roiId,iRun));
        ppiName.a = sprintf('R%i-%s*a',iRun,roiId);
        ppiName.b = sprintf('R%i-%s*b',iRun,roiId);
        %% Need to loop over runs and make all vois and then a ppi for each voi
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.spmmat = {spmMatFn};
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.adjust = 1;
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.session = iRun;
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.name = roiId;
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.roi{1}.mask.image = {roiFn};
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.roi{2}.mask.image = {epiMask};
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.roi{2}.mask.threshold = 0.5;
        spmBatch{batchIdxs(((iRun-1)*3 +1))}.spm.util.voi.expression = 'i1&i2';

        spmBatch{batchIdxs(((iRun-1)*3 +2))}.spm.stats.ppi.spmmat = {spmMatFn};
        spmBatch{batchIdxs(((iRun-1)*3 +2))}.spm.stats.ppi.type.ppi.voi = {voiFn};
        spmBatch{batchIdxs(((iRun-1)*3 +2))}.spm.stats.ppi.type.ppi.u = [1 1 1];
        spmBatch{batchIdxs(((iRun-1)*3 +2))}.spm.stats.ppi.name = ppiName.a ;
        spmBatch{batchIdxs(((iRun-1)*3 +2))}.spm.stats.ppi.disp = 0;

        spmBatch{batchIdxs(((iRun-1)*3 +3))}.spm.stats.ppi.spmmat = {spmMatFn};
        spmBatch{batchIdxs(((iRun-1)*3 +3))}.spm.stats.ppi.type.ppi.voi = {voiFn};
        spmBatch{batchIdxs(((iRun-1)*3 +3))}.spm.stats.ppi.type.ppi.u = [2 1 1]; %the first colum indexes SPM.Sess.U(i)
        spmBatch{batchIdxs(((iRun-1)*3 +3))}.spm.stats.ppi.name = ppiName.b ;
        spmBatch{batchIdxs(((iRun-1)*3 +3))}.spm.stats.ppi.disp = 0;
    end
    spm_jobman('initcfg');
    spm_jobman('run',spmBatch);
    toMove = [dirs.A02,filesep,'PPI_*-',roiId,'*'];
    movefile(toMove,dirs.PPI);
end
return
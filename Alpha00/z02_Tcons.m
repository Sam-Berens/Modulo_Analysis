function [] = z02_Tcons(G)
subjectIds = getSubjectIds(G);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
for iSubject=1:numel(subjectIds)

    cId = subjectIds{iSubject};
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Alpha00 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha00'];

    % Loop through the stimuli to compute map-wide t-statistics.
    % These contrasts sums scaled beta across runs
    for iStim = 0:5
        stimId = ['i',num2str(iStim)];
        spmMat = [dirs.Alpha00 ,filesep,stimId,filesep,'SPM.mat'];
        spmJob{1}.spm.stats.con.spmmat = {spmMat};
        spmJob{1}.spm.stats.con.consess{1}.tcon.name = stimId;

        % Weights is the unpadded weight vector.
        % It relies on the design matrix being organised such that the
        % stimulus of intrst is always encoded by the first column.
        spmJob{1}.spm.stats.con.consess{1}.tcon.weights = 1;
        spmJob{1}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
        spmJob{1}.spm.stats.con.delete = 1;
        spm_jobman('initcfg');
        spm_jobman('run',spmJob);
    end
end
return

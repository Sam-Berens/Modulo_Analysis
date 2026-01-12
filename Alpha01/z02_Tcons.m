function [] = z02_Tcons(G)
subjectIds = getSubjectIds(G);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
for iSubject = 1:numel(subjectIds)

    cId = char(subjectIds(iSubject));
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];

    % Loop through each position and stimulus to compute the t-stats
    pos = 'ab';
    for iPos = 1:2
        for iStim = 0:5
            stimId = sprintf('%s%i',pos(iPos),iStim);
            spmMat = [dirs.Alpha01,filesep,stimId,filesep,'SPM.mat'];
            spmJob{1}.spm.stats.con.spmmat = {spmMat};
            spmJob{1}.spm.stats.con.consess{1}.tcon.name = stimId;
            
            % Weights is the unpadded weight vector.
            % It relies on the design matrix being organised such that the
            % stimulus of intrest is always encoded by the first column.
            spmJob{1}.spm.stats.con.consess{1}.tcon.weights = 1;
            spmJob{1}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
            spmJob{1}.spm.stats.con.delete = 1;
            spm_jobman('initcfg');
            spm_jobman('run',spmJob);
        end
    end
end
return
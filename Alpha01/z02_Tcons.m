function [] = z02_Tcons(G)
subjectIds = getSubjectIds(G);
subjectIds = cellstr(subjectIds);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
for iSubject=1:numel(subjectIds)

    cId = subjectIds{iSubject};
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];

    %loop through the regs of intrst 0:5 and calculate their
    % tstatistics (in a contrast which scales the beta for each run
    pos = 'ab';
    for iStim=0:5
        for iPos=1:2
            stimId = sprintf('%s%i',pos(iPos),iStim);
            spmMat = [dirs.Alpha01,filesep,stimId,filesep,'SPM.mat'];
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
end
return

function [] = z02_Tcons(G)
subjectIds = getSubjectIds(G);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
for iSubject=1:numel(subjectIds)

    cId = subjectIds{iSubject};
    dirs.Subject = [dirs.Data,filesep,cId];
    dirs.Alpha01 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha01'];

    %loop through the regs of intrst 0:5 and calculate their
    % tstatistics (in a contrast which scales the beta for each run
    pos = 'ab';
    for ii=0:5
        for iPos=1:2
            rgName = sprintf('%s%i',pos(iPos),ii);
            spmMat = [dirs.Alpha01,filesep,rgName,filesep,'SPM.mat'];
            spmJob{1}.spm.stats.con.spmmat = {spmMat};
            spmJob{1}.spm.stats.con.consess{1}.tcon.name = 'rgName';
            %this is the unpadded weight vector, it relies on the beta
            %structure being organised such that the reg of intrst is always 1st
            spmJob{1}.spm.stats.con.consess{1}.tcon.weights = 1;
            spmJob{1}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
            spmJob{1}.spm.stats.con.delete = 1;
            spm_jobman('initcfg');
            spm_jobman('run',spmJob);
        end
    end
end
return

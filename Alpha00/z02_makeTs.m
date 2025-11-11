function [] = z02_makeTs(G)
subjectIds = getSubjectIds(G);
dirs.Data = ['..',filesep,'..',filesep,'Data'];
for ss=1:numel(subjectIds)

    cId = subjectIds{iSubject};
    dirs.Subject = [dirs.data,filesep,cId];
    dirs.Alpha00 = [dirs.Subject,filesep,'Analysis',filesep,'Alpha00'];

    %loop through the regs of intrst 0:5 and calculate their
    % tstatistics (in a contrast which scales the beta for each run
    for ii=0:5
        rgName = ['i',num2str(ii)];
        spmMat = [dirs.Alpha00 ,filesep,rgName,'SPM.mat'];
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
return

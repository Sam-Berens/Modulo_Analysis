function [] = singleDicom2Nifti(SubjectId, Filepath)

    cd(['Y:\Modulo\Data\',SubjectId])
    subjectDir = [pwd,filesep,'DICOM'];
    
    matlabbatch{1}.spm.util.import.dicom.data = {[subjectDir,Filepath]};
    matlabbatch{1}.spm.util.import.dicom.root = 'series';
    matlabbatch{1}.spm.util.import.dicom.outdir = {subjectDir};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.meta = 1;
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
    
    % Call the jobman and retun
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);

return
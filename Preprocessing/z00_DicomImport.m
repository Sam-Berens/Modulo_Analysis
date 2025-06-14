function [] = z00_DicomImport(SubjectId)

scriptsDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);

%% Setting the directory paths
workingDir = pwd;
subjectDir = [pwd,filesep,SubjectId];
dicomDir = [subjectDir,filesep,'DICOM'];

%% Get a cell array of all DICOMS
cd(dicomDir);
sessionName = dir();
sessionName = sessionName(3).name;
cd(sessionName);
fileList = dir();
fileList = fileList(3:end,1);
fileList = cellfun(@(s1,s2) [s1,filesep,s2],...
    {fileList.folder}',{fileList.name}','UniformOutput',false);
cd(workingDir);

%% Spec the job
spmBatch{1}.spm.util.import.dicom.data = fileList;
spmBatch{1}.spm.util.import.dicom.root = 'series';
spmBatch{1}.spm.util.import.dicom.outdir = {subjectDir};
spmBatch{1}.spm.util.import.dicom.protfilter = '.*';
spmBatch{1}.spm.util.import.dicom.convopts.format = 'nii';
spmBatch{1}.spm.util.import.dicom.convopts.meta = 1;
spmBatch{1}.spm.util.import.dicom.convopts.icedims = 0;

%% Call the jobman and return
spm_jobman('initcfg');
spm_jobman('run',spmBatch);

%% CD back to scripts
cd(scriptsDir);

return
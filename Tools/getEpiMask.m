function [epiMask] = getEpiMask(subjectId) 
if iscategorical(subjectId)
    subjectId = char(subjectId);
end
dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.Subject = [dirs.Data,filesep,subjectId];
dirs.EPI = [dirs.Subject,filesep,'EPI'];
epiMask.name = sprintf('%s%s_%s_epiMask00.nii',dirs.EPI,filesep,subjectId);
epiMask.V = spm_vol(epiMask.name);
epiMask.M = spm_read_vols(epiMask.V);
epiMask.idx = find(epiMask.M > 0.5);
epiMask.size = size(epiMask.M);
return
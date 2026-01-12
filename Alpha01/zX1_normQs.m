function [] = zX1_normQs()
dirs.Data = '../../Data';
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
for iSubject=1:nSubs
    subjectId = char(subjectIds(iSubject));
    trgtFld = fullfile(dirs.Data,subjectId,'Analysis','Alpha01',...
    'Mdl02');
    trgNames = {'q.nii','error.nii','b1.nii','b0.nii'};
    imagePaths = cellfun(@(x) fullfile(trgtFld,x), trgNames,...
        'UniformOutput',false);
    norm2MNI(G,subjectId,imagePaths');
end
return

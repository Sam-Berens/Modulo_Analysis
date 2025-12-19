function [] = zX1_normQs()
dirs.Data = '../../Data';
G = 'G1';
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
for iSubject=1:nSubs
    subjectId = subjectIds(iSubject);
    norm2MNI(G,subjectId);
end
return

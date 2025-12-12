function [nBackC] = getNback(G)
TaskIO = getScanTaskIO(G);
subjectIds = categorical(getSubjectIds(G));
nSubs = numel(subjectIds);
nBackC = nan(nSubs,1);
for iSubject = 1:nSubs
    sId = subjectIds(iSubject);
    sTaskIO = TaskIO(TaskIO.SubjectId==sId,:);
    s = (sTaskIO.TrialType=='1Back') | (sTaskIO.TrialType=='2Back');
    sNback = sTaskIO(s,:);
    nBackC(iSubject,1) = sum([sNback.correct]);
end
return

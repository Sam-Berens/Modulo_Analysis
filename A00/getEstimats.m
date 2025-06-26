function [Estimates] = getEstimats()
dirLits = dir('../../Data');
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectIds = categorical(subjectIds);
Estimates = [];
for iSubject = 1:numel(subjectIds)
    sId = subjectIds(iSubject);
    X = load(sprintf('..%s..%sData%s%s%sBehavioural%sA00.mat',...
        filesep,filesep,filesep,char(sId),filesep,filesep));
    Estimates = [Estimates;X.Estimates];
end
return
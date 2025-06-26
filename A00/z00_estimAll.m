function [] = z00_estimAll()
dirLits = dir('../../Data');
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectIds = categorical(subjectIds);
for iSubject = 1:numel(subjectIds)
    sId = subjectIds(iSubject);
    [~,lsRes] = system(sprintf(...
        'ls -l ..%s..%sData%s%s%sBehavioural%s | grep A00.mat',...
        filesep,filesep,filesep,char(sId),filesep,filesep));
    if contains(lsRes,'A00.mat')
        continue
    end
    t0 = table(...
        repmat(sId,4,1),...
        'VariableNames',{'SubjectId'});
    try
        [t1,Data] = estimFreq_spth01(char(sId),false);
        Estimates = [t0,t1];
        Data.RT = [];
        Data.FieldIdx_R = [];
        save(sprintf('..%s..%sData%s%s%sBehavioural%sA00.mat',...
            filesep,filesep,filesep,char(sId),filesep,filesep),...
            "-fromstruct",struct('Estimates',Estimates,'Data',Data));
        disp(iSubject);
    catch
        fprintf('Estimation failure for %s;%c',char(sId),10)
    end
end
return
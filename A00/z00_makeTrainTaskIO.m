function [] = z00_makeTrainTaskIO()
dirLits = dir('../../Data');
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectIds = categorical(subjectIds);
fh = waitbar(0,'Saving data...');
for iSubject = 1:numel(subjectIds)
    sId = subjectIds(iSubject);
    [~] = getTrainTaskIO(char(sId));
    waitbar(iSubject/numel(subjectIds),fh);
end
close(fh);
return
function [T] = getTotalTrainDuration(G)
T = get_pNonc(G);
T.TotalTrainDuration = nan(size(T,1),1);
dirs.Data = fullfile('..','..','Data');
for iSubject = 1:size(T,1)
    dirs.Subject = fullfile(dirs.Data,char(T.subjectId(iSubject)));
    dirs.Behavioural = fullfile(dirs.Subject,'Behavioural');
    TrainTaskIO = load(fullfile(dirs.Behavioural,'TrainTaskIO.mat'));
    T.TotalTrainDuration(iSubject) = ...
        sum(cellfun(@(v)(v(end)./1000)+1,TrainTaskIO.TaskIO.RT))/(60^2);
end
return
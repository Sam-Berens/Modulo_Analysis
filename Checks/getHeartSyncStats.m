function [HeartSync] = getHeartSyncStats(G)
TaskIO = getScanTaskIO(G);
SubjectId = unique(TaskIO.SubjectId);
statsToDo = {'phiKeyp1','phiRespo'};
figure;
for iStat = 1:numel(statsToDo)
    n = nan(size(SubjectId));
    c = nan(size(SubjectId));
    mag = nan(size(SubjectId));
    arg = nan(size(SubjectId));
    z = nan(size(SubjectId));
    pVal = nan(size(SubjectId));
    subplot(1,2,iStat);
    hold on;
    for iSubject = 1:numel(SubjectId)
        T = TaskIO(TaskIO.SubjectId==SubjectId(iSubject),:);
        theta = T.(statsToDo{iStat});
        theta = theta(~isnan(theta));
        n(iSubject) = numel(theta);
        [pVal(iSubject), z(iSubject)] = circ_rtest(theta);
        c(iSubject) = mean(exp(1i*theta));
        mag(iSubject) = abs(c(iSubject));
        arg(iSubject) = angle(c(iSubject));
        plot([0,real(c(iSubject))],[0,imag(c(iSubject))]);
    end
    xlim([-0.2,0.2]);
    ylim([-0.2,0.2]);
    axis square;
    title(statsToDo{iStat});
    HeartSync.(statsToDo{iStat}) = table(SubjectId,n,c,mag,arg,z,pVal);
    fprintf('%s%c',statsToDo{iStat},10)
    disp(HeartSync.(statsToDo{iStat}));
    disp('');
end
return
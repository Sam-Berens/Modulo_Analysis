function [HeartSync] = getHeartSyncStats(G)
TaskIO = getScanTaskIO(G);
subjectId = unique(TaskIO.SubjectId);
statsToDo = {'phiKeyp1','phiRespo'};
figure;
for iStat = 1:numel(statsToDo)
    n = nan(size(subjectId));
    c = nan(size(subjectId));
    mag = nan(size(subjectId));
    arg = nan(size(subjectId));
    z = nan(size(subjectId));
    pVal = nan(size(subjectId));
    subplot(1,2,iStat);
    hold on;
    for iSubject = 1:numel(subjectId)
        T = TaskIO(TaskIO.SubjectId==subjectId(iSubject),:);
        theta = T.(statsToDo{iStat});
        theta = theta(~isnan(theta));
        n(iSubject) = numel(theta);
        if n(iSubject) > 0
            [pVal(iSubject), z(iSubject)] = circ_rtest(theta);
            c(iSubject) = mean(exp(1i*theta));
            mag(iSubject) = abs(c(iSubject));
            arg(iSubject) = angle(c(iSubject));
            plot([0,real(c(iSubject))],[0,imag(c(iSubject))]);
        end
    end
    xlim([-0.2,0.2]);
    ylim([-0.2,0.2]);
    axis square;
    title(statsToDo{iStat});
    HeartSync.(statsToDo{iStat}) = table(subjectId,n,c,mag,arg,z,pVal);
    fprintf('%s%c',statsToDo{iStat},10)
    disp(HeartSync.(statsToDo{iStat}));
    disp('');
end
return
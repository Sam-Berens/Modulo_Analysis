function [dt03] = getDataTable03()
fname = 'DataTable03.mat';
if exist(fname,'file')
    strct = load(fname);
    dt03 = strct.dt03;
    return
end

%set a seed
rng(1729);

G = 'G1';
taskIO = getScanTaskIO(G);
taskIO.Properties.VariableNames{1} = 'subjectId';
heartO = getHeartSyncStats(G);
dt03 = join(taskIO,heartO.phiKeyp1);
%filter out the null trials
dt03(isnan(dt03.iTrial),:) = [];
%Wipe the null category from memory
dt03.TrialType = categorical(cellstr(dt03.TrialType));

%for any missed responses, we are going to inject a correct response 1/6
%times
nNans = sum(isnan(dt03.correct),1);
dt03{isnan(dt03.correct),'correct'} = binornd(1,1/6,nNans,1);

save(fname,'dt03');
return;

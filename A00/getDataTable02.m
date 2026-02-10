function [DataTable02] = getDataTable02()
G = 'G1';
subjectIds = getSubjectIds(G);
DataTable00 = load('DataTable00.mat');
DataTable00 = DataTable00.DataTable00;
DataTable00 = DataTable00(ismember(DataTable00.subjectId,subjectIds),:);
HeartStats = getHeartSyncStats(G);
DataTable02 = join(DataTable00,HeartStats.phiKeyp1);
return
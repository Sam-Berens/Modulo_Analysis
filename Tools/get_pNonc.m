function [pNonc] = get_pNonc(G)
% pNonc is a vector within a table encoding the summed probability of a
% correct response for all non-commutable pairs at the end of training, as
% estimated by taking the mean posterior values from A00/vonMises.

model = 'vonMises';
dirs.Data = ['..',filesep,'..',filesep,'Data'];
dirs.Group = [dirs.Data,filesep,'_Group'];
dirs.A00 = [dirs.Group,filesep,'A00'];
dirs.vonMises = [dirs.A00,filesep,model];

% Load the data
Output = load([dirs.vonMises,filesep,'Output.mat'],'Profic');
Profic = Output.Profic;
Profic = Profic.EAP;

pairIds = [8,13,17,21,28,32]; % These are zero-ordered!
colIdx = pairIds + 2; % One to account for subjectId, one for 1-ordering
Profic = Profic(:,[1,colIdx]);

% Create a table
subjectId = Profic.subjectId;
pNonc = sum(table2array(Profic(:,2:end)),2);
pNonc = table(subjectId,pNonc);

% Sort/select the rows to ensure we are consistent with G
subjectId = getSubjectIds(G);
ii = arrayfun(@(c)find(c==pNonc.subjectId),subjectId);
pNonc = pNonc(ii,:);
return
function [mcp_noncom] = getP_noncom(G)

%𝑝_𝑛𝑜𝑛𝑐𝑜𝑚 is a mean-centred predictor encoding the summed probability ...
% of a correct response for all non-commutable pairs at the end of training,..
% as estimated by A00 using the Expected A Posteriori.

%for now defualt:
model = 'vonMises';
subjectIds = getSubjectIds(G);
[~,order] = sort(subjectIds);

cd ../../Data/_Group/A00
cd(model);
strct = load('Output.mat');

profs = strct.Profic;
%load in the Expected A Posteriori for each PairId, as estimed by A00
allProfs = profs.EAP;

nSubs = numel(subjectIds);
p_noncom = nan(nSubs,1);
for iSubject = 1:nSubs
    cId = subjectIds(iSubject);
    cProfs = allProfs(allProfs.subjectId == cId,:);
    %cut off the subjectId column so indexing is clearer
    cProfs = cProfs(1,2:end);
    %select out the columns of the table which are non com
    nonComIds = [8,13,17,21,28,32]; %we zero ordered them elsewhere
    %columns are not zero ordered...
    nonComIds = nonComIds +1;
    nonComPrfs = cProfs{1,nonComIds};
    p_noncom(iSubject,1) = sum(nonComPrfs);
end



%mean centre the predictor
mcp_noncom = p_noncom- mean(p_noncom);
%make sure its sorted by subjectId
mcp_noncom = mcp_noncom(order,1);

%return to home of this script
cd  ../../../../Scripts/Alpha01/

return

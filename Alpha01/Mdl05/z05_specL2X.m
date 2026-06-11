function [] = z05_specL2X(G)

% Set dirs
dirs.Data = fullfile('..','..','..','Data');
dirs.Output = fullfile(dirs.Data,...
    '_Group',...
    G,...
    'Analysis',...
    'Alpha01',...
    'Mdl05',...
    'GLM0');

% Get subject Ids
subjectId = getSubjectIds(G);
nSubjects = numel(subjectId);

% Get subject effects
subjectFxs = pivot(...
    table(subjectId),'columns','subjectId','rows','subjectId');

% Get zPnonc
pNonc = get_pNonc(G);
pNonc.zPnonc = zscore(pNonc.pNonc);
zPnonc = removevars(pNonc,"pNonc");

% Make X
X = join(subjectFxs,zPnonc);
X = repelem(X,2,1);
X.colocation = repmat([-1;1],nSubjects,1);
X = removevars(X,'subjectId');
X = movevars(X,'zPnonc','after','colocation');
X.('zPnonc:colocation') = X.zPnonc.*X.colocation;
X.zPnonc = []; % Remove zPnonc as it is perfectly correlated with subject

% Save the result
R = X{:,:};
names = X.Properties.VariableNames;
outputFn = fullfile(dirs.Output,'X.mat');
if ~exist(dirs.Output,"dir")
    mkdir(dirs.Output);
end
save(outputFn,'R','names');

return

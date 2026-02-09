function [] = z05_make2ndLvlDesign()
% Cd out into Alpha01
wd = pwd;
cd ..;
dirs.Data = '../../Data';
G = 'G1';
subjectId = getSubjectIds(G);
nSubs = numel(subjectId);
tWide = pivot(table(subjectId),'columns','subjectId', ...
    'rows','subjectId');%doing like this to keep subjects in the right order
t = get_pNonc(G);
%t = sortrows(t,'subjectId');
zPnonc = zscore(t.pNonc);
pNonc = [t,table(zPnonc)];
zPnonc = removevars(pNonc,"pNonc");
tWide = join(tWide,zPnonc);
%tWide = pivot(pNonc,'columns','subjectId', 'rows','zPnonc');
tWide = repelem(tWide,2,1);
colocation = nan(nSubs*2,1);
for ii = 1:(nSubs*2)
    % iSub = ceil((ii/2));
    if rem(ii,2)
        colocation(ii) = -1; %note that this will be in the first row for each subject
    else
        colocation(ii) = +1;
    end
end
X = [tWide,table(colocation)];
X = movevars(X,'zPnonc','after','colocation');
X = removevars(X,'subjectId'); %this is just taking off the 'long' version not the 'wide'

names = X.Properties.VariableNames;
names = [names,{'zPnonc:colocation'}];
R = X{:,:};

%add interaction column;
xTermStruct = getXtermStruct(names);
R = addXterms(R,xTermStruct);
fldr = fullfile(dirs.Data,'_Group',G,'Analysis','Alpha01',...
    'Mdl04');
fname = fullfile(fldr,'X.mat');

if ~exist(fldr,"dir")
    mkdir(fldr)
end 
save(fname,'R','names');

%Cd back in 
cd(wd);

return

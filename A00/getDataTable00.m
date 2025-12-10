function [DataTable00] = getDataTable00()
%Constructs table of subject, pairId, qType, Chi-term, learnt predictor,etc

if exist('DataTable00.mat','file')
    DataTable00 = load('DataTable00.mat');
    DataTable00 = DataTable00.DataTable00;
    return
end

dirs.model = fullfile('..','..','Data','_Group','A00','vonMises');

stanOut = load([dirs.model,filesep,'Output.mat']);
% Expected A Posteriori for each questions a and b parameter in the model
Estims = stanOut.Estims.EAP(:,[1,18:end]);
varNames = Estims.Properties.VariableNames;

s = contains(varNames,'a_');
s(1) = true; % Keep subjectID
A = Estims(:,s);
A = Wide2Long(A,'a');

s = contains(varNames,'b_');
s(1) = true; % Keep subjectID
B = Estims(:,s);
B = Wide2Long(B,'b');

s = contains(varNames,'c_');
s(1) = true; % Keep subjectID
C = Estims(:,s);
C = Wide2Long(C,'c');

% Pr(r>0.5 @ maxX)
Prg05 = stanOut.Profic.prg05;
Prg05 = Wide2Long(Prg05,'prg05');


DataTable00 = join(A,B);
DataTable00 = join(DataTable00,C);
DataTable00 = join(DataTable00,Prg05);

% Make learnt
DataTable00.learnt = double(DataTable00.prg05 > 0.5);
%everything beyond this is post-meeting addtions
[names, codes] = getQtype(Prg05.pairId);
DataTable00.pairType = names; %I think technically the model function can cope with using this instead of the binary columns but check
DataTable00.sup = codes(:,1); 
DataTable00.zero = codes(:,2); 
DataTable00.commute = codes(:,3); 
DataTable00.noncom = codes(:,4); 
%Make Chi term 𝜒
sSup = DataTable00.pairType=="Sup";
DataTable00.chi = nan(size(DataTable00,1),1);
DataTable00.chi(sSup,1) = ...
    DataTable00.a(sSup) + DataTable00.c(sSup);
DataTable00.chi(~sSup,1) = DataTable00.a(~sSup);
save('DataTable00.mat',"DataTable00");
return



function [pairTypeName,typeBool] = getQtype(pairId)
% normalize input to a column vector
if iscategorical(pairId)
    %we're doing it this way because turning a categorical array ...
    %straight to a double as different behaviour ...
    %depending on if input is scalar or array char(pairId);
    pairId = char(pairId);
    pairId = str2num(pairId);
end
pairId = pairId(:);
if any(pairId>35)
    error('invalid pairId')
end
allIdx = (0:35)';
zer = [(0:5)'; (6:6:30)'];
unC = [3;6;10;23;26;31;34];
unC = unC(~ismember(unC,zer));
unN = [8,13,17,21,28,32]';
unN = unN(~ismember(unN,zer));
sup = allIdx(~ismember(allIdx,[zer;unC;unN]));
% make horizontal cell so cellfun returns cells corresponding to columns
typeDict = {sup, zer, unC, unN};
% each cellfun call returns an Nx1 logical; keep as cells then convert
typeBoolCell = cellfun(@(x) ismember(pairId, x), typeDict, 'UniformOutput', false);
typeBool = cell2mat(typeBoolCell);  % Nx4 logical matrix
% use max to get the column index (type) for each row
[~, pairCode] = max(typeBool, [], 2);
%make them into string names
typeNames = categorical({'Sup', 'Zer', 'UnC', 'UnN'});
pairType = arrayfun(@(x) typeNames(x),pairCode,'UniformOutput',false);
pairTypeName = cellfun(@(x) x, pairType);
return

function Long = Wide2Long(Wide,newVarName)
Wide.Properties.VariableNames(2:end) = ...
    arrayfun(@(ii){sprintf('%02d',ii)},0:35);
varsToStack = Wide.Properties.VariableNames(2:end);
Long = stack(Wide, varsToStack, ...
    'NewDataVariableName',newVarName, ...
    'IndexVariableName','pairId');
Long.pairId = categorical(Long.pairId);
return
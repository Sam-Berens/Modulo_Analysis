function [DataTable00] = getDataTable00()
%Constructs table of subject, pairId, qType, Chi-term, learnt predictor,etc

if exist('DataTable00.mat','file')
    DataTable00 = load('DataTable00.mat');
    DataTable00 = DataTable00.DataTable00;
    return
end

dirs.model = fullfile('..','..','Data','_Group','A00','vonMises');
stanOut = load(fullfile(dirs.model,'Output.mat'));

% MAP/EAP tables, trimmed
Estims = {stanOut.Estims.MAP(:,[1,18:end]), ...
    stanOut.Estims.EAP(:,[1,18:end])};
%extract maxX and preallocate 
MaxX = stanOut.MaxX;
DataTable00 = table();
measures = ["MAP","EAP"]; %types of estimates to store for each term 
params   = ["a","b","c"];  % prefixes to extract

% Pre-extract pNan (only used for c)
pNanRaw = stanOut.Estims.pNan;
pNanVarNames = pNanRaw.Properties.VariableNames;

% Convert pNan to long format once per measure later
cMask_pNan = contains(pNanVarNames,"c_");
cMask_pNan(1) = true;  % keep subjectID
pNanWide = pNanRaw(:,cMask_pNan);

for iM = 1:2
    Est = Estims{iM};
    measure = measures(iM);
    vNames = Est.Properties.VariableNames;
    Results = struct();
    % Handle a, b, c in same loop 
    for p = params
        mask = contains(vNames, p + "_");
        mask(1) = true;  % keep subjectID
        wideTable = Est(:,mask);

        % Convert wide to long
        longTable = Wide2Long(wideTable, measure + "_" + p);
        % post-processing for "c" using proportion of samples which were
        % either predicted to be complex of inf in A00 
        if p == "c"
            % join MaxX + pNan
            LT = join(longTable, MaxX);
            LT = join(LT, Wide2Long(pNanWide, "pNan"));
            % less negatively biased c estimate
            LT.ubc = LT.pNan .* 2.*LT.maxX + LT.(measure + "_c") .* (1-LT.pNan);
            % rename to avoid duplicate columns later
            LT.Properties.VariableNames{'pNan'} = char(measure + "c_pNan");
            LT.Properties.VariableNames{'ubc'} = char(measure + "ubc");
            longTable = LT;
        end
        Results.(p) = longTable;
    end

    % append running table
    if iM == 1
        DataTable00 = join(Results.a, Results.b);
        DataTable00 = join(DataTable00, Results.c);
        % Pr(r>0.5 @ maxX)
        Prg05 = stanOut.Profic.prg05;
        Prg05 = Wide2Long(Prg05,'prg05');
        DataTable00 = join(DataTable00,Prg05);
        % Make learnt
        DataTable00.learnt = double(DataTable00.prg05 > 0.5);
        [names, codes] = getQtype(Prg05.pairId);
        DataTable00.pairType = names;
        DataTable00.sup = codes(:,1);
        DataTable00.zero = codes(:,2);
        DataTable00.commute = codes(:,3);
        DataTable00.noncom = codes(:,4);
    else
        DataTable00 = join(DataTable00, Results.a);
        DataTable00 = join(DataTable00, Results.b);
        DataTable00 = join(DataTable00, Results.c);
    end
    %Make Chi term for either MAP or EAP
    sSup = DataTable00.pairType=="Sup";
    DataTable00.(measure+"chi") = nan(size(DataTable00,1),1);
    DataTable00.(measure+"chi")(sSup,1) = ...
    DataTable00.(measure + "_a")(sSup) + DataTable00.(measure + "_c")(sSup);
    DataTable00.(measure+"chi")(~sSup,1) = DataTable00.(measure + "_a")(~sSup);
end
%move around columns created by weird looping artifacts
DataTable00 = movevars(DataTable00, 'maxX', 'After', 'EAPchi');
DataTable00 = movevars(DataTable00, 'MAPchi', 'After', 'MAPubc');

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

function [DesignMat]=getDesignMat(G)
%finds the model-predicted n feedback events for r=0.5

dirs.home = pwd;
cd ../../Data/_Group/A00/vonMises
dirs.vM = pwd;
cd(dirs.home);

fn = [dirs.vM,filesep,'Output.mat'];
output = load(fn);
%Expected A Posteriori for each questions a and b parameter in the model
Estims = output.Estims.EAP;
Estims = sortrows(Estims,'subjectId'); 
%Pr(r>0.5 @ maxTrial) 
Prg05 = output.Profic.prg05;
Prg05 = sortrows(Prg05,'subjectId');
Prg05 = Prg05(:,2:end); %check that chopping off subject id works ok

subjectIds = categorical(getSubjectIds(G));
nSubs = numel(subjectIds);
nQs = 36;

%build columns of design matrix
Chi = nan((nSubs*nQs),1);
Learnt = nan((nSubs*nQs),1);
Zero = zeros((nSubs*nQs),1);
Commute = zeros((nSubs*nQs),1);
Noncom = zeros((nSubs*nQs),1);
PairId = nan((nSubs*nQs),1);
SubjectId = categorical((nSubs*nQs),1);

%set our specified r value
%r = pC2r(0.99);
r = pC2r(0.5);

for ii=1:numel(subjectIds)
    cId = subjectIds(ii,1);
    for jj=1:36
        qId = jj-1;
        qType = getQtype(qId);

        varNm = ['a_',int2str(jj)];
        a = Estims{ii,varNm}; 

        varNm = ['c_',int2str(jj)];
        c = Estims{1,varNm};

        cRow = ((ii-1)*36)+jj;
        %add current pair pre-table vector
        PairId(cRow,1) = jj; %should i make this zero orderered?
        %contruct Chi term 𝜒𝑗,𝑠={
        % 𝑎𝑗,𝑠 + 𝜏𝑗,𝑠,    
        if isequal(qType,1)
            chi = a + c;
         %𝑗∈ 𝑆𝑢𝑝   𝑎𝑗,s,
        else
            chi = a;
            %add current Qtype to pre-table vectors
            if isequal(qType,1)
                Zero(cRow,1)= 1;
            elseif isequal(qType,2)
                Commute(cRow,1)= 1;
            else
                Noncom(cRow,1)= 1;
            end
        end
        
        %this is wrong - we want the probability that r is greater than 0.5
        %so we can check that against 0.5, not whether the r value
        %predicted at the end of training was bigger than 0.5
        cPrg05 = Prg05{ii,jj};
        learnt = cPrg05 > 0.5;
        
        %add to Data and predictors to pre-table vectors
        Chi(cRow,1) = chi;
        Learnt(cRow,1) = learnt;
        SubjectId(cRow,1) = cId;
    end
end
PairId = categorical(PairId);

DesignMat = table(Chi,Learnt,Zero,Commute,Noncom,PairId,SubjectId);
return



function [qType] = getQtype(qId)

allIdx = (0:35)';
zer = [(0:5)';(6:6:30)'];
unC = [3;6;10;23;26;31;34];
unC = unC(~ismember(unC,zer));
unN = [8,13,17,21,28,32]';
unN = unN(~ismember(unN,zer));
sup = allIdx(~ismember(allIdx,[zer;unC;unN]));

typeDict = {sup;zer;unC;unN};
typeBool = cellfun(@(x) ismember(qId,x),typeDict);
qType = find(typeBool);

return


function [XMem] = getXmem(G)
%finds the model-predicted n feedback events for r=0.5

dirs.home = pwd;
cd ../../Data/_Group/A00/vonMises
dirs.vM = pwd;
cd(dirs.home);

fn = [dirs.vM,filesep,'Output.mat'];
output = load(fn);
estims = output.Estims.EAP;
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
nQs = 36;

%build column of desing matrix
XMem = nan((nSubs*nQs),1);

%set our specified r value
r = pC2r(0.99);

for ii=1:numel(subjectIds)
    cId = subjectIds{ii};
    maxX = getMaxX(cId);
    alpha.i1 = estims{ii,'alpha1_1'};
    alpha.i2 = estims{ii,'alpha1_2'};
    alpha.i3 = estims{ii,'alpha1_3'};
    alpha.i4 = estims{ii,'alpha1_4'};
    beta.i1 =  estims{ii,'beta1_1'};
    beta.i2 =  estims{ii,'beta1_2'};
    beta.i3 =  estims{ii,'beta1_3'};
    beta.i4 =  estims{ii,'beta1_4'};
    
    alpha.s1 = estims{ii,'alpha2_1'};
    alpha.s2 = estims{ii,'alpha2_2'};
    alpha.s3 = estims{ii,'alpha2_3'};
    alpha.s4 = estims{ii,'alpha2_4'};
    beta.s1 =  estims{ii,'beta2_1'};
    beta.s2 =  estims{ii,'beta2_2'};
    beta.s3 =  estims{ii,'beta2_3'};
    beta.s4 =  estims{ii,'beta2_4'};

    for jj=1:36
        qId = jj-1;
        qType = getQtype(qId);
        s_qType = ['s',int2str(qType)];
        i_qType = ['i',int2str(qType)];

        varNm = ['a_',int2str(jj)];
        alphaSlope = estims{ii,varNm};
        alphaSlope = alphaSlope + alpha.(s_qType);
        alphaInt = alpha.(i_qType);
        a = sigmoid((alphaInt + alphaSlope)) + maxX;   %TO DO - something about what i'm givning sigmoid function is wrong!    
        
        varNm = ['b_',int2str(jj)];
        betaSlope = estims{ii,varNm};
        betaSlope = betaSlope + beta.(sType);
        betaInt = beta.(i_qType);
        b = betaInt + betaSlope;

        xMem = log(exp((atanh(r)/b))-1) + a;
        cRow =((ii-1)*36)+jj;
        XMem(cRow,1) = xMem;
    end
end

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

function [maxX] = getMaxX(subjectId)
% Make paths structure
dir.Data = ['..',filesep,'..',filesep,'Data',filesep];
dir.Sub = [dir.Data,filesep,subjectId];
dir.Beh = [dir.Sub,filesep,'Behavioural'];

strct = load([dir.Beh,filesep,'TrainTaskIO.mat']);
taskIO = strct.TaskIO;
% x: Learning state predictor, Real[0,Inf]
x = taskIO.tSup;
maxX = max(x);
return
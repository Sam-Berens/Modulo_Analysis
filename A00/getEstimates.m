function [] = getEstimates(subjectIds)

% Turn off warnings
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
warning('off','MATLAB:print:ExcludesUIInFutureRelease');

%% Set the params to extract
param = [];
for ii = 1:4
    param = [param,string(sprintf('alpha1_%i',ii))];
    param = [param,string(sprintf('alpha2_%i',ii))];
    param = [param,string(sprintf('beta1_%i',ii))];
    param = [param,string(sprintf('beta2_%i',ii))];
end
for ii = 1:36
    param = [param,string(sprintf('a_%i',ii))];
    param = [param,string(sprintf('b_%i',ii))];
end

%% Make dirPaths structure
dirPaths = struct;
dirPaths.Data = ['..',filesep,'..',filesep,'Data',filesep];
dirPaths.IO = @(sId) [dirPaths.Data,sId,filesep,...
    'Analysis',filesep,'A00',filesep];

dirPaths.groupOut = [dirPaths.Data,filesep,'_Group',filesep,'A00'];
if ~exist(dirPaths.groupOut,'dir')
    mkdir(dirPaths.groupOut);
end
dirPaths.groupFigs = [dirPaths.groupOut,filesep,'Figures'];
if ~exist(dirPaths.groupFigs,'dir')
    mkdir(dirPaths.groupFigs);
end

%% Get subjectIds
if nargin == 0
    dirLits = dir(dirPaths.Data);
    subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
end
subjectIds = categorical(subjectIds);

%% Subject loop
for iSubject = 1:numel(subjectIds)
    sId = char(subjectIds(iSubject));
    inFn = [dirPaths.IO(sId),'InputData.json'];
    outList = dir([dirPaths.IO(sId),'Output_*.csv']);

    % Set the domain
    inData = jsondecode(fileread(inFn));
    maxX = inData.maxX;

    % Get posterior samples
    P = readCSVs(sId,outList);

    % Plot the postrior samples
    plotStanFit(P,dirPaths.groupFigs);

    % Extract proficiency
    Q = expractProfic(sId,P,maxX);

    % Extract estimates
    E = extractEstims(sId,P,param,maxX);

    % Get the diagnostics
    [R,D] = getDiagnostics(sId,P,param);

    % Construct the big tables
    if iSubject == 1
        PostSamps = P;
        Profic = Q;
        Estims = E;
        Rhat = R;
        Diverg = D;
        Pfns = fieldnames(Profic);
        Efns = fieldnames(Estims);
    else
        PostSamps = [PostSamps;P];
        Rhat = [Rhat;R];
        Diverg = [Diverg;D];
        for iF = 1:numel(Pfns)
            Profic.(Pfns{iF}) = [Profic.(Pfns{iF});Q.(Pfns{iF})];
        end
        for iF = 1:numel(Efns)
            Estims.(Efns{iF}) = [Estims.(Efns{iF});E.(Efns{iF})];
        end
    end
end

%% Print Rhat
plotRhat(Rhat,dirPaths.groupFigs);

%% Turn warning back on
warning('on','MATLAB:table:ModifiedAndSavedVarnames');
warning('on','MATLAB:print:ExcludesUIInFutureRelease');

%% Package diagnostics and save
Diagno.Rhat = Rhat;
Diagno.Diverg = Diverg;

save([dirPaths.groupOut,filesep,'Output.mat'],...
    'PostSamps',...
    'Profic',...
    'Estims',...
    'Diagno');

return

function [StanOut] = readCSVs(subjectId,outList)
for iChain = 1:size(outList,1)
    filePath = [outList(iChain).folder,filesep,outList(iChain).name];
    lines = readlines(filePath);
    s = false(size(lines,1),1);
    for ii = 1:numel(s)
        if ~isempty(lines{ii})
            c0 = lines{ii}(1);
            s(ii) = ~strcmp(c0,'#');
        end
    end
    lines = lines(s);
    tempFn = [tempname(pwd),'.csv'];
    writelines(lines,tempFn);
    T = readtable(tempFn,'VariableNamingRule','modify');
    delete(tempFn);
    chain = ones(size(T,1),1).*iChain;
    T = [table(chain),T]; %#ok<*AGROW>
    if iChain == 1
        StanOut = T;
    else
        StanOut = [StanOut;T];
    end
end
S = table(repmat(categorical({subjectId}),size(StanOut,1),1),...
    'VariableNames',{'subjectId'});
StanOut = [S,StanOut];
return

function [Profic] = expractProfic(subjectId,P,maxX)
theta = (0:5).*(pi/3);
[a,b] =  ind2sub([6,6],1:(6^2));
names = cellfun(@(i1,i2)sprintf('%i_%i',i1,i2),...
    num2cell(a-1),num2cell(b-1),'UniformOutput',false);

P0 = nan(size(P,1),6^2);
for ii = 1:(6^2)
    a = P.(['a_',num2str(ii)]);
    b = P.(['b_',num2str(ii)]);
    r = tanh(b.*log(1+exp(maxX-a)));
    k = r2k(r);
    pmf = exp(cos(theta).*k);
    pmf = pmf./sum(pmf,2);
    P0(:,ii) = pmf(:,1);
end

X = mat2cell(P0,size(P,1),ones(1,36));
Profic.EAP = cellfun(@getEAP,X);
Profic.MAP = cellfun(@getMAP,X);
Profic.pg0 = cellfun(@getPG0,X);
[Profic.low,Profic.med,Profic.upp] = cellfun(@getLMU,X);

fns = fieldnames(Profic);
for iF = 1:numel(fns)
    T = array2table(Profic.(fns{iF}),'VariableNames',names);
    Profic.(fns{iF}) = [...
        table(categorical({subjectId}),'VariableNames',{'subjectId'}),...
        T];
end
return

function [EAP] = getEAP(x)
EAP = mean(x);
return

function [x_map, details] = getMAP(x)
x = x(:);
N = numel(x);

% Choose a small window width h (rule-of-thumb from ranks; scale-free).
% This makes the method stable even if data pile up at 0/1.
h = max(1/(2*sqrt(N)), 1/N);  % Conservative small window

% Interior candidate via KDE on interior points only.
xin = x(x>0 & x<1);
x_mode_in = NaN;
if ~isempty(xin)
    % Keep KDE simple since we don't rely on its absolute scale
    [f,xi] = ksdensity(xin,'Support',[0,1],...
        'BoundaryCorrection','reflection');
    [~,i] = max(f);
    x_mode_in = xi(i);
end

% Small-ball probabilities
m0 = mean(x <= h);
m1 = mean((1-x) <= h);
mmid = 0;
if ~isnan(x_mode_in)
    mmid = mean(abs(x - x_mode_in) <= h);
end

% Pick the location with largest small-ball probability
[~,which] = max([m0, mmid, m1]);
if which == 1
    x_map = 0;
elseif which==3
    x_map = 1;
else
    x_map = x_mode_in;
end

if nargout > 1
    details = struct('h',h,...
        'm0',m0,...
        'mmid',mmid,...
        'm1',m1,...
        'x_mode_interior',x_mode_in);
end
return

function [pg0] = getPG0(x)
pg0 = mean(x>(1/6));
return

function [low,med,upp] = getLMU(x)
[f,xi] = ecdf(x);
[~,ii] = min((f-0.025).^2);
low = xi(ii);
[~,ii] = min((f-0.5).^2);
med = xi(ii);
[~,ii] = min((f-0.975).^2);
upp = xi(ii);
return

function [k] = r2k(r)
% r2k converts the inverse link function output values to a von Mises
% concentration parameter, kappa.
%
% Inputs:
%    r - A vector of Softplus tanh transformed predictor values.
%
% Outputs:
%    k - A vector of concentration parameters corresponding to r.
%
% The conversion uses piecewise definitions based on abs(r).
%
k = nan(size(r));
for ii = 1:1:numel(r)
    cr = abs(r(ii));
    if (cr < 0.85) && (cr >= 0.53)
        ck = -0.4 + (1.39*cr) + (0.43/(1-cr));
    elseif cr < 0.53
        ck = (2*cr) + (cr^3) + (5*(cr^5)/6);
    else
        ck = ((cr^3) - (4*(cr^2)) + (3*cr))^-1;
    end
    k(ii) = ck*sign(r(ii));
end
k(k>700) = 700;
k(k<-700) = -700;
return

function [Estims] = extractEstims(subjectId,StanOut,param,maxX)
Estims.EAP = nan(1,numel(param));
Estims.MAP = nan(1,numel(param));
Estims.pg0 = nan(1,numel(param));
Estims.low = nan(1,numel(param));
Estims.med = nan(1,numel(param));
Estims.upp = nan(1,numel(param));

isSd = find(contains(param,'2_'));
isA = find(contains(param,'a_'));
for iParam = 1:numel(param)
    x = StanOut.(param(iParam));

    % Expected value
    Estims.EAP(iParam) = mean(x);

    % Mode
    if ismember(iParam,isSd)
        [f,xi] = ksdensity(x,'Support','nonnegative',...
            'BoundaryCorrection','reflection');
        [~,imax] = max(f);
        Estims.MAP(iParam) = xi(imax);
    elseif ismember(iParam,isA)
        [f,xi] = ksdensity(x,'Support',[0,maxX],...
            'BoundaryCorrection','reflection');
        [~,imax] = max(f);
        Estims.MAP(iParam) = xi(imax);
    else
        [f,xi] = ksdensity(x);
        [~,imax] = max(f);
        Estims.MAP(iParam) = xi(imax);
    end

    % Probability greater than zero
    Estims.pg0(iParam) = mean(x>0);

    % 2.5%, Median, 97.5%
    [f,xi] = ecdf(x);
    [~,ii] = min((f-0.025).^2);
    Estims.low(iParam) = xi(ii);
    [~,ii] = min((f-0.5).^2);
    Estims.med(iParam) = xi(ii);
    [~,ii] = min((f-0.975).^2);
    Estims.upp(iParam) = xi(ii);
end
fns = fieldnames(Estims);
for iF = 1:numel(fns)
    T = array2table(Estims.(fns{iF}),'VariableNames',param);
    Estims.(fns{iF}) = [...
        table(categorical({subjectId}),'VariableNames',{'subjectId'}),...
        T];
end
return

function [Rhat,Diverg] = getDiagnostics(subjectId,StanOut,param)

% Rhat
SumStats = join(...
    groupsummary(StanOut,'chain','mean',param),...
    groupsummary(StanOut,'chain','var',param));
N = SumStats.GroupCount(1);
M = SumStats.Variables;
idxMeans = cellfun(...
    @(s)find(strcmp(s,SumStats.Properties.VariableNames)),...
    cellfun(@(s)sprintf('mean_%s',s),cellstr(param),...
    'UniformOutput',false));
idxVars = cellfun(...
    @(s)find(strcmp(s,SumStats.Properties.VariableNames)),...
    cellfun(@(s)sprintf('var_%s',s),cellstr(param),...
    'UniformOutput',false));
B = N.*var(M(:,idxMeans));
W = mean(M(:,idxVars),1);
V =((N-1)/N).*W + ((1/N).*B);
Rhat = sqrt(V./W);
Rhat = [table(categorical({subjectId}),'VariableNames',{'subjectId'}),...
    array2table(Rhat,'VariableNames',param)];

% Divergences
SumStats = groupsummary(StanOut,'chain','sum',"divergent__");
chain = SumStats.chain;
count = SumStats.sum_divergent__;
total = sum(count);
Diverg = array2table([count',total],...
    'VariableNames',[...
    cellfun(@(ii)sprintf('chain_%i',ii),num2cell(chain),...
    'UniformOutput',false);{'total'}]);
Diverg = [...
    table(categorical({subjectId}),'VariableNames',{'subjectId'}),...
    Diverg];

return

function [] = plotStanFit(StanOut,figPath)
subjectId = char(StanOut.subjectId(1));
pointGroup = StanOut.chain;
figPath = [figPath,filesep,subjectId];
if ~exist(figPath,'dir')
    mkdir(figPath);
end

% Types
types = {'ZeroPlus';'Supervised';'UnsCom';'UnsNon'};
for iT = 1:numel(types)
    Ax = ['alpha1_',num2str(iT)]; % Location
    Ay = ['alpha2_',num2str(iT)]; % Scale
    Bx = ['beta1_',num2str(iT)]; % Location
    By = ['beta2_',num2str(iT)]; % Scale

    fh = figure('Units','normalized','OuterPosition',[0,0,1,1]);
    u1 = uipanel(fh,'position',[.0,0,.5,1],'Title','alpha');
    u2 = uipanel(fh,'position',[.5,0,.5,1],'Title','beta');
    scatterhist(StanOut.(Ax),StanOut.(Ay),...
        'Group',pointGroup,...
        'Kernel','on','Parent',u1);
    xlabel('Mean');
    ylabel('Std');
    sgtitle([types{iT},': Alpha']);
    scatterhist(StanOut.(Bx),StanOut.(By),...
        'Group',pointGroup,...
        'Kernel','on','Parent',u2);
    xlabel('Mean');
    ylabel('Std');
    sgtitle([types{iT},': Beta']);

    print(fh,sprintf('%s%sG%i_%s.png',...
        figPath,filesep,iT,types{iT}),'-dpng');
    close(fh);
end

% Pairs
for iP = 1:36
    fh = figure('Units','normalized','OuterPosition',[0,0,1,1]);
    [iy,ix] = ind2sub([6,6],iP);
    scatterhist(...
        StanOut.(['a_',num2str(iP)]),...
        StanOut.(['b_',num2str(iP)]),...
        'Group',pointGroup,...
        'Kernel','on');
    xlabel('a');
    ylabel('b');
    sgtitle(sprintf('%i+%i',(ix-1),(iy-1)));

    print(fh,sprintf('%s%sP%i_%i+%i.png',...
        figPath,filesep,iP,(ix-1),(iy-1)),'-dpng');
    close(fh);
end
return

function [] = plotRhat(Rhat,savePath)
fh = figure('Units','normalized','OuterPosition',[0,0,1,1]);
R = Rhat(:,2:end).Variables;
imagesc(R);
colorbar;
xticks(1:size(R,2));
xtickangle(90);
xticklabels(strrep(Rhat.Properties.VariableNames(2:end),'_','.'));
yticks(1:size(Rhat,1));
yticklabels(Rhat.subjectId);
print(fh,sprintf('%s%sRhat.png',...
    savePath,filesep),'-dpng');
close(fh);
return
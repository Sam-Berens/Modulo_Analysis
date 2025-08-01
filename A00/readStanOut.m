function [Rhat,StanOut] = readStanOut()


dirList = dir('cc*.csv');
for iFile = 1:size(dirList,1)
    lines = readlines(dirList(iFile).name);
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
    T = readtable(tempFn);
    delete(tempFn);
    chain = ones(size(T,1),1).*iFile;
    T = [table(chain),T]; %#ok<*AGROW>
    if iFile == 1
        StanOut = T;
    else
        StanOut = [StanOut;T];
    end
end

%%
p = figure;
u1 = uipanel(p,'position',[0,.5,.5,.5],'Title','ZeroPlus');
u2 = uipanel(p,'position',[.5,.5,.5,.5],'Title','Supervised');
u3 = uipanel(p,'position',[0,0,.5,.5],'Title','UnsCom');
u4 = uipanel(p,'position',[.5,0,.5,.5],'Title','UnsNon');
scatterhist(StanOut.alpha1_1,StanOut.beta1_1,'Group',StanOut.chain,'Kernel','on','Parent',u1);
scatterhist(StanOut.alpha1_2,StanOut.beta1_2,'Group',StanOut.chain,'Kernel','on','Parent',u2);
scatterhist(StanOut.alpha1_3,StanOut.beta1_3,'Group',StanOut.chain,'Kernel','on','Parent',u3);
scatterhist(StanOut.alpha1_4,StanOut.beta1_4,'Group',StanOut.chain,'Kernel','on','Parent',u4);

%%
pp = figure('Units','normalized','OuterPosition',[0,0,1,1]);
for iPlot = 1:36
    [iy,ix] = ind2sub([6,6],iPlot);
    pos = [(ix-1)/6,(-1/6)+1-((iy-1)/6),1/6,1/6];
    uu = uipanel(pp,'position',pos,'Title',sprintf('%i+%i',(ix-1),(iy-1)));
    scatterhist(StanOut.(['a_',num2str(iPlot)]),StanOut.(['b_',num2str(iPlot)]),'Group',StanOut.chain,'Kernel','on','Parent',uu);
end

%% Diagnotsics
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
diagnost = join(...
    groupsummary(StanOut,'chain','mean',param),...
    groupsummary(StanOut,'chain','var',param));
N = diagnost.GroupCount(1);
temp = diagnost.Variables;
idxMeans = cellfun(...
    @(s)find(strcmp(s,diagnost.Properties.VariableNames)),...
    cellfun(@(s)sprintf('mean_%s',s),cellstr(param),...
    'UniformOutput',false));
idxVars = cellfun(...
    @(s)find(strcmp(s,diagnost.Properties.VariableNames)),...
    cellfun(@(s)sprintf('var_%s',s),cellstr(param),...
    'UniformOutput',false));
B = N.*var(temp(:,idxMeans));
W = mean(temp(:,idxVars),1);
V =((N-1)/N).*W + ((1/N).*B);
Rhat = sqrt(V./W)';
Rhat = table(param',Rhat);

%%
figure;
histogram(StanOut.divergent__);

return
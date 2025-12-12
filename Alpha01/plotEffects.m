function [] = plotEffects(roiName)
%get normal design matrix
DataTable00 = getDataTable00;

%select out apt roi
dX = DataTable00(...
    DataTable00.roiName==['l',roiName] | ...
    DataTable00.roiName==['r',roiName], :);
formula = 'zTemplate ~ coLocation*mCpNonc*hemisphere + (1|subjectId)';
mdl = fitlme(dX,formula);

meanPnonc = mean(dX(1:31,7));
% hlfdomain = 6- meanPnonc{1,1};
% domain = linspace(hlfdomain,6)';
% domain = domain - meanPnonc{1,1};
domain = linspace(-2,4); 

%return the fixed effect predictors coding table
dX = designMatrix(mdl, 'Fixed');

terms = mdl.CoefficientNames;
terms = terms';

idx1 =  find(contains(terms,'coLocation'));
left = nan(1,numel(terms));
left(1,1) = 1; %intercept
left(1,idx1) = 0; %we're taking the middle of the coLocation effect (because its a bipolar predictor)
idx2 = find(strcmp(terms,'hemisphere'));
left(1,idx2) = 1;
idx3 = find(contains(terms,'mCpNonc'));
idx4 = idx3(~ismember(idx3,idx1)); %find fixed effect and hemi interaction only
left(1,idx4) = NaN; %this is so we can put in the domain later

panel1 = repmat(left,[numel(domain),1]);
panel1(isnan(panel1)) = [domain,domain]; %there are 2 nan columns which is why we have to do this

panel2 = panel1.*-1 ; %this makes the coding represent left hemi instead
panel2(:,1) = 1 ;%make the intercept positive again

%now we make hemi effects 0 and coLocation 1 (just opposite of above)
coLoc = nan(1,numel(terms));
coLoc(1,1) = 1; %intercept
idx5 = find(contains(terms,'hemisphere'));
coLoc(1,idx5) = 0 ;
idx6 = find(strcmp(terms,'coLocation'));
coLoc(1,idx6) = 1; %coLoc = aa/bb comparisons
idx7 = idx3(~ismember(idx3,idx5));%this is to get the pnonc and pnonc:coloc
coLoc(1,idx7) = NaN;
panel3 = repmat(coLoc,[numel(domain),1]);
panel3(isnan(panel3)) = [domain,domain];

panel4 =  panel3.*-1 ; %this makes the coding represent left hemi instead
panel4(:,1) = 1 ;%make the intercept positive again

panels = {panel1;panel2;panel3;panel4};
pNonc = figure;
titles = {'leftHemi';'rightHemi';'coLocated';'heteroLocated'};
for ii=1:4
in = num2cell(panels{ii,:},1); 
[~,est,low,upp] = getXEUL(mdl,in{:});
%add the subplot
subplot(2,2,ii);
plot(domain,est,'DisplayName','est');
hold on;
plot(domain,low,'DisplayName','low');
plot(domain,upp,'DisplayName','upp');
legend;
xlabel('mCpNonc');
ylabel('zTemplate');
title(titles{ii,1});
hold off;
end 

ax = findall(gcf,'type','axes');   % get all axes handles
set(ax, 'YLim', [-0.1 0.3]);   % set common limits

%interactions would be:
%you want a design mat where the interaction is on and off but there are
%like 3 different ways the interaction can be off?



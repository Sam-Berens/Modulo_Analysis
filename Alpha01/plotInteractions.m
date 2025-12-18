function [] = plotInteractions(mdlOb,dataTb)

meanPnonc = mean(dataTb(1:31,7));
% hlfdomain = 6- meanPnonc{1,1};
% domain = linspace(hlfdomain,6)';
% domain = domain - meanPnonc{1,1};
domain = linspace(-2,4); 


terms = mdlOb.CoefficientNames;
terms = terms';

%panel 1 is going to be coLoc =1 hemi=1
idxs1 =  find(contains(terms,'coLocation'));
panel1 = nan(1,numel(terms));
panel1(1,1) = 1; %intercept
panel1(1,idxs1) = 1; %coLoc aa/bb comparisons
idxs2 = find(contains(terms,'hemisphere'));
panel1(1,idxs2) = 1;
%Tag anything with pnonc as nan so we can put in the domain 
idx3 = find(contains(terms,'mCpNonc'));
panel1(1,idx3) = NaN; 

panel1 = repmat(panel1,[numel(domain),1]);
%there are 4 nan columns which is why we have to do this
panel1(isnan(panel1)) = repmat(domain,[1,4]); 
%panel 2 is going to be coLoc =1 hemi=-1
panel2 = panel1;
panel2(:,idxs2) = panel2(:,idxs2).*-1 ; %turn all hemi containing terms negative

%panel 3 is going to be coLoc=-1 hemi=1
panel3 = nan(1,numel(terms));
panel3(1,1) = 1; %intercept
panel3(1,idxs2) = 1 ;% hemi 
panel3(1,idxs1) = -1; %coLoc = ab comparisons
panel3(1,idx3) = NaN; %now all pnonc involved terms are tagged for domain
panel3 = repmat(panel3,[numel(domain),1]);
panel3(isnan(panel3)) =  repmat(domain,[1,4]); 
%sets any interactions containing coLoc:mcPnonc to be negative again)
idx4 = find(contains(terms,'coLocation:mCpNonc'));
panel3(:,idx4) = panel3(:,idx4).*-1;

%panel 4 is going to be coLoc=-1 hemi=-1 
panel4 =  panel1; 
panel4(:,[2,4,5,7]) =  panel4(:,[2,4,5,7]).*-1 ;

panels = {panel1;panel2;panel3;panel4};
pNonc = figure;
% titles = {'leftHemi';'rightHemi';'coLocated';'heteroLocated'};
for ii=1:4
in = num2cell(panels{ii,:},1); 
[~,est,low,upp] = getXEUL(mdlOb,in{:});
%add the subplot
subplot(2,2,ii);
plot(domain,est,'DisplayName','est');
hold on;
plot(domain,low,'DisplayName','low');
plot(domain,upp,'DisplayName','upp');
legend;
xlabel('mCpNonc');
ylabel('zTemplate');
% title(titles{ii,1});
hold off;
end 
subplot(2,2,1) 
title('Left');
subplot(2,2,2)
title('Right');
% Row labels
annotation('textbox',[0.05 0.45 0.05 0.2], ...
    'String','AA/BB comparison', ...
    'FontWeight','bold',...
    'FitBoxToText','on',...
    'LineStyle','none',...
    'Rotation',90)

annotation('textbox',[0.05 0.01 0.05 0.2], ...
    'String','AB comparison', ...
    'FontWeight','bold',...
    'FitBoxToText','on',...
    'LineStyle','none',...
    'Rotation',90)


ax = findall(gcf,'type','axes');   % get all axes handles
set(ax, 'YLim', [-0.1 0.3]);   % set common limits
return